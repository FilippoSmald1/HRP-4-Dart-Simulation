#include "MPCSolver.hpp"
#include "utils.cpp"
#include "qpOASES/qpOASES.hpp"
#include  <stdlib.h>
#include <iostream>
#include <fstream>
#include  <Eigen/Cholesky>

using namespace mpcSolver;

MPCSolver::MPCSolver(double mpcTimeStep, double controlTimeStep, double predictionTime, Eigen::Vector3d initialComPosition,
                     double comTargetHeight, double singleSupportDuration, double doubleSupportDuration, double thetaMax,
					 double footConstraintSquareWidth, double deltaXMax, double deltaYIn, double deltaYOut, double measuredComWeight, double measuredZmpWeight, bool _activate_timing_adaptation ){

	// Set up parameters
	this->mpcTimeStep = mpcTimeStep;
	this->controlTimeStep = controlTimeStep;

	this->footConstraintSquareWidth = footConstraintSquareWidth;
	this->deltaXMax = deltaXMax;
	this->deltaYIn = deltaYIn;
	this->deltaYOut = deltaYOut;
	this->thetaMax = thetaMax;

        
   
    this->activate_timing_adaptation = _activate_timing_adaptation;
    this->comTargetHeight = comTargetHeight;
    
    this->omega=sqrt(9.81/(comTargetHeight));


    this->singleSupportDuration=singleSupportDuration;
    this->doubleSupportDuration=doubleSupportDuration;

    this->measuredComWeight_x = measuredComWeight;
    this->measuredZmpWeight = measuredZmpWeight;

    measuredComWeight_y = measuredComWeight_x;

    N = round(predictionTime/mpcTimeStep);
    S = round(singleSupportDuration/mpcTimeStep);
    D = round(doubleSupportDuration/mpcTimeStep);
    M = ceil(N/(S+D));
    mpcIter = 0;
    controlIter = 0;

    // Timing Manager Initialization
    // The matrix is composed by threee columns with the following roles
    // ->> timing shrink, gait phase*, timing buffer
    // * single support = 1, doublesupport = 0
    Timing_Manager = Eigen::MatrixXd::Zero(6,3);
    Timing_Manager << singleSupportDuration, 1, singleSupportDuration, doubleSupportDuration, 0, doubleSupportDuration, singleSupportDuration, 1, singleSupportDuration, doubleSupportDuration, 0, doubleSupportDuration, 0, 1, 0, 0, 0, 0;

    // Matrices for cost function
    costFunctionH = Eigen::MatrixXd::Zero(2*(N+M),2*(N+M));
    costFunctionF = Eigen::VectorXd::Zero(2*(N+M));

    // Matrices for stability constraint
    Aeq = Eigen::MatrixXd::Zero(2,(N*2)+(M*2));
    beq = Eigen::VectorXd::Zero(2);

    // Matrices for ZMP constraint
    AZmp = Eigen::MatrixXd::Zero(2*N,2*(N+M));
    bZmpMax = Eigen::VectorXd::Zero(2*N);
    bZmpMin = Eigen::VectorXd::Zero(2*N);

    // Matrices for feasibility constraint
    AFootsteps = Eigen::MatrixXd::Zero(2*M,2*(M+N));
    bFootstepsMax = Eigen::VectorXd::Zero(2*M);
    bFootstepsMin = Eigen::VectorXd::Zero(2*M);

    // Matrices for swing foot constraint
    ASwingFoot = Eigen::MatrixXd::Zero(4,(N*2)+(M*2));
    bSwingFoot = Eigen::VectorXd::Zero(4);

    // Matrices for all constraints stacked
    AConstraint = Eigen::MatrixXd::Zero(2*(N+M)+2,2*(N+M));
    bConstraintMin = Eigen::VectorXd::Zero(2*(N+M)+2);
    bConstraintMax = Eigen::VectorXd::Zero(2*(N+M)+2);

    // Matrices for ZMP prediction
    p =  Eigen::VectorXd::Ones(N);
    P =  Eigen::MatrixXd::Ones(N,N)*mpcTimeStep;

    for(int i=0; i<N;++i){
    	for(int j=0;j<N;++j){
            if (j>i) P(i,j)=0;
        }
    }

    // Matrices for CoM velocity prediction and LIP state update
    Vu = Eigen::MatrixXd::Zero(N,N);
    Vs = Eigen::MatrixXd::Zero(N,3);

    double ch = cosh(omega*mpcTimeStep);
    double sh = sinh(omega*mpcTimeStep);

    Eigen::MatrixXd A_upd = Eigen::MatrixXd::Zero(3,3);
    Eigen::VectorXd B_upd = Eigen::VectorXd::Zero(3);
    A_upd<<ch,sh/omega,1-ch,omega*sh,ch,-omega*sh,0,0,1;
    B_upd<<mpcTimeStep-sh/omega,1-ch,mpcTimeStep;

    Eigen::RowVectorXd Vu_newline(N);
    Eigen::RowVectorXd Vs_newline(3);
    Eigen::RowVectorXd A_midLine(3);

    A_midLine<<omega*sh,ch,-omega*sh;

    for(int i=1;i<=N;++i) {
        Vu_newline.setZero();
        Vs_newline.setZero();
        Vu_newline(i-1) = 1-ch;
        if(i>1) {
            for(int j=0;j<=i-2;++j) {
                Vu_newline.segment(j,1) = A_midLine*(matrixPower(A_upd,(i-j-2)))*B_upd;
            }
        }
        Vs_newline = A_midLine*(matrixPower(A_upd,(i-1)));
        Vu.row(i-1) = Vu_newline;
        Vs.row(i-1) = Vs_newline;
    }

    // Vector that contains predicted rotations
    predictedOrientations = Eigen::VectorXd::Zero(M+1);

    // Initialize CoM, ZMP and predicted footstep
    comPos = initialComPosition;
    comPos(2) = comTargetHeight;
    comVel = Eigen::Vector3d::Zero(3);
    zmpPos = Eigen::Vector3d(comPos(0),comPos(1),0.0);

    predictedFootstep = Eigen::Vector4d::Zero();

    // Initialize footstep counter
    footstepCounter=0;

    // A storage matrix
    predictedZmp = Eigen::MatrixXd::Zero(3,N);


}

void MPCSolver::solve(Eigen::Vector3d measuredComPos, Eigen::Vector3d measuredComVel, Eigen::Vector3d measuredComAcc,
					  Eigen::Affine3d swingFootTransform, bool supportFoot, double simulationTime, double vRefX, double vRefY, double omegaRef, bool widgetReference){

    // Save iteration parameters

if (!widgetReference) {
    this->vRefX = vRefX;
    this->vRefY = vRefY;
    this->omegaRef = omegaRef;
}

    this->supportFoot = supportFoot;
    this->simulationTime = simulationTime;

    //std::cout<<Timing_Manager<<std::endl;




    // If new footstep is starting, change reference frame
    if(supportFootHasChanged()) changeReferenceFrame(swingFootTransform);
    //if(Timing_Manager(5,0) == 1) changeReferenceFrame(swingFootTransform);
  
    // Adjust the state based on measures
 
    if (InitCom == 0){

    measuredComWeight_x = 1;
    measuredComWeight_y = 1;
    measuredComWeight_v_x = 0; 
    measuredComWeight_v_y = 0;
    measuredZmpWeight = 0.0;

    InitCom = InitCom + 1;

}else{

    if(footstepCounter > 3){ 

    // DCM check 
    if (pow((measuredComPos(0) +  measuredComVel(0)/omega)-(comPos(0)+comVel(0)/omega),2)>0.0005 && trig_x ){ //>0.001
        measuredComWeight_x = 0;    
        measuredComWeight_v_x = 1;
        trig_x = false;
        trig_y = false;
    }else{ 

    measuredComWeight_v_x = 0;
    measuredComWeight_x = 0;
}
    if (pow((measuredComPos(1) +  measuredComVel(1)/omega)-(comPos(1)+comVel(1)/omega),2)>0.001 && trig_y ) {//>0.001
         measuredComWeight_y = 0;  
         measuredComWeight_v_y = 1;
         trig_x = false;
         trig_y = false;
 }   else{ 

    measuredComWeight_v_y = 0;
    measuredComWeight_y = 0;
}

    measuredZmpWeight = 0.0;

   
}else{ 

    measuredComWeight_x = 0.0;
    measuredComWeight_y = 0.0;
    measuredComWeight_v_x = 0; 
    measuredComWeight_v_y = 0;
    measuredZmpWeight = 0.0;

}

}
     


/*
        // Uncomment to understand the performance of the kinematic controller:
        // measCom and MPCCom (and measVel and comVel) should be almost the same, when no perturbations or reference velocity changes occur! 
        std::cout <<"measCom " << measuredComPos <<std::endl;
        std::cout <<"MPCCom " << comPos <<std::endl;
        std::cout <<"measVel " << measuredComVel <<std::endl;
        std::cout <<"comVel " << comVel <<std::endl;
/**/

    comPos(0) = (1-measuredComWeight_x)*comPos(0) + measuredComWeight_x*measuredComPos(0);
    comPos(1) = (1-measuredComWeight_y)*comPos(1) + measuredComWeight_y*measuredComPos(1);

    comVel(0) = (1-measuredComWeight_v_x)*comVel(0) + measuredComWeight_v_x*measuredComVel(0);
    comVel(1) = (1-measuredComWeight_v_y)*comVel(1) + measuredComWeight_v_y*measuredComVel(1);

    Eigen::Vector3d measuredZmpPos = measuredComPos - measuredComAcc/(omega*omega);
    zmpPos(0) = (1-measuredZmpWeight)*zmpPos(0) + measuredZmpWeight*measuredZmpPos(0);
    zmpPos(1) = (1-measuredZmpWeight)*zmpPos(1) + measuredZmpWeight*measuredZmpPos(1);



           

  


    // Compute footsteps orientations
    computeOrientations();

    // Compute the cost function
    genCostFunction();

    // Compute the matrices for the constraints
    genStabilityConstraint();
    genBalanceConstraint();
    genFeasibilityConstraint();
    genSwingFootConstraint(swingFootTransform);

    // Stack the matrices for the constraints
    // If both feet are on the floor also add the swing foot constraint

    if (false) { //(footstepCounter == 0 || mpcIter >= 0)
    	int nConstraints = Aeq.rows() + AFootsteps.rows() + AZmp.rows() + ASwingFoot.rows();
    	AConstraint.resize(nConstraints, 2*(N+M));
    	bConstraintMin.resize(nConstraints);
    	bConstraintMax.resize(nConstraints);

    	AConstraint 	 << Aeq, AFootsteps, AZmp, ASwingFoot;
    	bConstraintMin  << beq, bFootstepsMin, bZmpMin, bSwingFoot;
    	bConstraintMax << beq, bFootstepsMax, bZmpMax, bSwingFoot;
    } else {
    	int nConstraints = Aeq.rows() + AFootsteps.rows() + AZmp.rows();
    	AConstraint.resize(nConstraints, 2*(N+M));
    	bConstraintMin.resize(nConstraints);
    	bConstraintMax.resize(nConstraints);

    	AConstraint 	 << Aeq, AFootsteps, AZmp;
    	bConstraintMin  << beq, bFootstepsMin, bZmpMin;
    	bConstraintMax << beq, bFootstepsMax, bZmpMax;
    }

    // Solve QP
    Eigen::VectorXd decisionVariables = solveQP();
    
    // Split the QP solution in ZMP dot and footsteps
    Eigen::VectorXd zDotOptimalX(N);
    Eigen::VectorXd zDotOptimalY(N);
    Eigen::VectorXd footstepsOptimalX(M);
    Eigen::VectorXd footstepsOptimalY(M);

    zDotOptimalX = (decisionVariables.head(N));
    zDotOptimalY = (decisionVariables.segment(N+M,N));
    footstepsOptimalX = decisionVariables.segment(N,M);
    footstepsOptimalY = decisionVariables.segment(2*N+M,M);

    // Save some stuff for plotting
    predictedZmp.row(0) = (P*zDotOptimalX + p*zmpPos(0)).transpose();
    predictedZmp.row(1) = (P*zDotOptimalY + p*zmpPos(1)).transpose();
    predictedZmp.row(2) = Eigen::RowVectorXd::Zero(N);

    // Update the state based on the result of the QP
    Eigen::Vector3d nextStateX = updateState(zDotOptimalX(0),0,controlTimeStep);
    Eigen::Vector3d nextStateY = updateState(zDotOptimalY(0),1,controlTimeStep);

    comPos << nextStateX(0),nextStateY(0),comTargetHeight;
    comVel << nextStateX(1),nextStateY(1),0.0;
    zmpPos << nextStateX(2),nextStateY(2),0.0;
    predictedFootstep << footstepsOptimalX(0),footstepsOptimalY(0),0.0,predictedOrientations(1);


    ++controlIter;


    mpcIter = floor(controlIter*controlTimeStep/mpcTimeStep);


  // Timing manager update
     if (Timing_Manager(0,0)-controlTimeStep >= 0.0) {

         Timing_Manager(0,0) = Timing_Manager(0,0) - controlTimeStep;
         Timing_Manager(4,0) = Timing_Manager(4,0) + controlTimeStep;
         Timing_Manager(4,2) = Timing_Manager(0,2);
         Timing_Manager(4,1) = Timing_Manager(0,1);
         Timing_Manager(5,0) = 0;

     }else{
         
        trig_y = true;
        trig_x = true;

     if (Timing_Manager(4,1) == doublesupport) Timing_Manager(5,0) = 1;

        Eigen::MatrixXd line = Eigen::MatrixXd::Zero(1,3);
        line << 0, 1, 0;
        Timing_Manager.block(0,0,1,3) = Timing_Manager.block(1,0,1,3);
        Timing_Manager.block(1,0,1,3) = Timing_Manager.block(2,0,1,3);
        Timing_Manager.block(2,0,1,3) = Timing_Manager.block(3,0,1,3);
        Timing_Manager.block(3,0,1,3) = Timing_Manager.block(4,0,1,3);
        Timing_Manager(3,0) = Timing_Manager(3,0) + controlTimeStep; 
        Timing_Manager.block(4,0,1,3) = line;
     }

    // Updating countdown of adapted step
    if (CountDown != -100){

    if (CountDown -1 == 2){
    // Reset timing manager matrix, with footstep change (since the adapted step has been performed)
    Timing_Manager << singleSupportDuration, 1, singleSupportDuration, doubleSupportDuration, 0, doubleSupportDuration, singleSupportDuration, 1, singleSupportDuration, doubleSupportDuration, 0, doubleSupportDuration, 0, 1, 0, 1, 0, 0; 
    CountDown = -100;
    }else CountDown = CountDown-1;

    }



if(activate_timing_adaptation){
    if(Timing_Manager(5,0)==1){
        trig_y = true;
        trig_x = true;
    	controlIter = 0;
    	mpcIter = 0;
        footstepCounter++;
    }
}else{
    if(mpcIter>=S+D){ 
        trig_y = true;
        trig_x = true;
    	controlIter = 0;
    	mpcIter = 0;
        footstepCounter++;
    }

}

  


    std::cout << "Iteration " << controlIter << " Footstep " << footstepCounter << std::endl;
}


void MPCSolver::changeReferenceFrame(Eigen::Affine3d swingFootTransform) {

	Eigen::Vector3d swingFootPos = swingFootTransform.translation();

	// apply rotation before translation because this transform
	// is expressed in the new support foot frame
	comPos = swingFootTransform.rotation()*comPos;
	comVel = swingFootTransform.rotation()*comVel;
	zmpPos = swingFootTransform.rotation()*zmpPos;


	comPos(0)-= -swingFootPos(0);
	comPos(1)-= -swingFootPos(1);
	zmpPos(0)-= -swingFootPos(0);
	zmpPos(1)-= -swingFootPos(1);
}

void MPCSolver::genStabilityConstraint() {
	Eigen::VectorXd b(N);

	for(int i=0;i<N;++i){
		b(i) = pow(exp(-omega*mpcTimeStep),i);
	}

	Aeq.block(0,0,1,N)   = ((1/omega)*((1-exp(-omega*mpcTimeStep))/(1+pow(exp(-omega*mpcTimeStep),N))))*b.transpose();
	Aeq.block(1,N+M,1,N) = ((1/omega)*((1-exp(-omega*mpcTimeStep))/(1+pow(exp(-omega*mpcTimeStep),N))))*b.transpose();

	beq<<comPos(0) + comVel(0)/omega - zmpPos(0), comPos(1) + comVel(1)/omega - zmpPos(1);
}

void MPCSolver::genSwingFootConstraint(Eigen::Affine3d swingFootTransform) {

	Eigen::Vector3d swingFootPos = swingFootTransform.translation();

	ASwingFoot(0,N) = 1;
	ASwingFoot(1,2*N+M) = 1;
	ASwingFoot(2,N+1) = 1;
	ASwingFoot(3,2*N+M+1) = 1;

	double step = 0.15;
	if (footstepCounter == 0) step = 0;

	int sign = 1;
	if (footstepCounter%2 == 0) sign = -1;

	bSwingFoot(0) = 0*step;//predictedFootstep(0);//swingFootPos(0);
	bSwingFoot(1) = sign*0.15;//predictedFootstep(1);//swingFootPos(1);
	bSwingFoot(2) = 0*2*step;
	bSwingFoot(3) = 0;
}

void MPCSolver::genCostFunction() {
	//qVx = 0;
	//qVy = 0;

	Eigen::MatrixXd Cc = Eigen::MatrixXd::Zero(N,M);
	Eigen::VectorXd Ccf = Eigen::VectorXd::Ones(S+D);

	for(int i=0; i<M-1; ++i){
		Cc.block(S+D-mpcIter+(i*(S+D)),i,S+D,1) = Ccf;
	}

	Cc.block(S+D-mpcIter+((M-1)*(S+D)),M-1,mpcIter,1) = Ccf.block(0,0,mpcIter,1);

	costFunctionH.block(0,0,N,N) = qZd*Eigen::MatrixXd::Identity(N,N) + qVx*Vu.transpose()*Vu + qZ*P.transpose()*P;
	costFunctionH.block(N+M,N+M,N,N) = qZd*Eigen::MatrixXd::Identity(N,N) + qVy*Vu.transpose()*Vu + qZ*P.transpose()*P;

	costFunctionH.block(0,N,N,M) = -qZ*P.transpose()*Cc;
	costFunctionH.block(N,0,M,N) = -qZ*Cc.transpose()*P;
	costFunctionH.block(N,N,M,M) = qZ*Cc.transpose()*Cc;

	costFunctionH.block(N+M,2*N+M,N,M) = -qZ*P.transpose()*Cc;
	costFunctionH.block(2*N+M,N+M,M,N) = -qZ*Cc.transpose()*P;
	costFunctionH.block(2*N+M,2*N+M,M,M) = qZ*Cc.transpose()*Cc;

	Eigen::VectorXd vArcX = Eigen::VectorXd::Zero(N);
	Eigen::VectorXd vArcY = Eigen::VectorXd::Zero(N);

	Eigen::VectorXd costFunctionF1 = Eigen::VectorXd::Zero(N+M);
	Eigen::VectorXd costFunctionF2 = Eigen::VectorXd::Zero(N+M);

	for(int i=0;i<N;++i){
		vArcX(i) = vRefX*cos(i*omegaRef*mpcTimeStep) - vRefY*sin(i*omegaRef*mpcTimeStep);
		vArcY(i) = vRefX*sin(i*omegaRef*mpcTimeStep) + vRefY*cos(i*omegaRef*mpcTimeStep);
	}

    Eigen::Vector3d stateX = Eigen::Vector3d(comPos(0),comVel(0),zmpPos(0));
    Eigen::Vector3d stateY = Eigen::Vector3d(comPos(1),comVel(1),zmpPos(1));

    costFunctionF1.block(0,0,N,1) = (qVx*Vu.transpose())*((Vs*stateX)-vArcX) + qZ*P.transpose()*p*zmpPos(0);
    costFunctionF2.block(0,0,N,1) = (qVy*Vu.transpose())*((Vs*stateY)-vArcY) + qZ*P.transpose()*p*zmpPos(1);

    costFunctionF1.block(N,0,M,1) = -qZ*Cc.transpose()*p*zmpPos(0);
    costFunctionF2.block(N,0,M,1) = -qZ*Cc.transpose()*p*zmpPos(1);

    costFunctionF<<costFunctionF1,costFunctionF2;
}

void MPCSolver::computeOrientations() {


	// Difference matrix (theta_j - theta_{j-1})
	Eigen::MatrixXd differenceMatrix = Eigen::MatrixXd::Identity(M,M);
	for (int i=0; i<M-1; ++i) {
		differenceMatrix(i+1,i) = -1;
	}

	// Rotation cost function
	Eigen::MatrixXd footstepsH = (differenceMatrix.transpose()*differenceMatrix);
	Eigen::VectorXd footstepsF = -omegaRef*(singleSupportDuration+doubleSupportDuration)*differenceMatrix.transpose()*Eigen::VectorXd::Ones(M);

	// Constraint on maximum variation of orientation
	Eigen::MatrixXd AOrientations = differenceMatrix;
	Eigen::VectorXd bOrientationsMin = -thetaMax*Eigen::VectorXd::Ones(M);
	Eigen::VectorXd bOrientationsMax =  thetaMax*Eigen::VectorXd::Ones(M);

	// Constraint on first step (no change in orientation)
	Eigen::MatrixXd AFirstOrientation = Eigen::MatrixXd::Zero(1,M);
	AFirstOrientation(0,0) = 1;
	Eigen::VectorXd bFirstOrientation = Eigen::VectorXd::Zero(M);

	// Some QP Options
	qpOASES::Options optionsOrientations;
	optionsOrientations.setToMPC();
	optionsOrientations.printLevel=qpOASES::PL_NONE;
	qpOASES::int_t nWSR = 300;

	qpOASES::QProblem qpRotations(M,M);
	qpRotations.setOptions(optionsOrientations);

	qpOASES::real_t thetaOpt[M];

	qpOASES::real_t footstepsHqpOASES[M*M];
	qpOASES::real_t footstepsFqpOASES[M];

	for(int i=0;i<M;++i){
		for(int j=0;j<M;++j){
			footstepsHqpOASES[i*M+j] = footstepsH(i,j);
		}
		footstepsFqpOASES[i] = footstepsF(i);
	}

	qpOASES::real_t AOrientationsQP[M*M];
	qpOASES::real_t bOrientationsMinQP[M];
	qpOASES::real_t bOrientationsMaxQP[M];

	for(int i=0;i<M;++i){
		for(int j=0;j<M;++j){
			AOrientationsQP[i*M+j] = AOrientations(i,j);
		}
		bOrientationsMinQP[i] = bOrientationsMin(i);
		bOrientationsMaxQP[i] = bOrientationsMax(i);
	}

	qpRotations.init(footstepsHqpOASES,footstepsFqpOASES,
			AOrientationsQP,0,0,bOrientationsMinQP,bOrientationsMaxQP,nWSR,NULL,NULL,NULL,NULL,NULL,NULL);

	qpRotations.getPrimalSolution(thetaOpt);

	predictedOrientations(0)=0;

	if(footstepCounter==0){
		for(int i=1;i<M;++i){
			predictedOrientations(1) = 0;
			predictedOrientations(i+1) = thetaOpt[i-1];
		}
	}
	else{
		for(int i=1;i<M+1;++i){
			predictedOrientations(i) = thetaOpt[i-1];
		}
	}
}

void MPCSolver::genBalanceConstraint(){
	AZmp.setZero();

        // Some local matrices and parameter initialization
	Eigen::MatrixXd Icf = Eigen::MatrixXd::Zero(S+D,S+D);
	Eigen::MatrixXd Ic = Eigen::MatrixXd::Zero(N,N);
	Eigen::MatrixXd Cc = Eigen::MatrixXd::Zero(N,M);
	Eigen::VectorXd Ccf = Eigen::VectorXd::Ones(S+D);
	Eigen::MatrixXd rCosZmp = Eigen::MatrixXd::Zero(N,N);
	Eigen::MatrixXd rSinZmp = Eigen::MatrixXd::Zero(N,N);
	Eigen::MatrixXd zmpRotationMatrix(2*N,2*N);

	for(int i=0;i<S;++i){
		Icf(i,i)=1;
	}

        int S_1,D_1,S_2,D_2,S_3,D_3;

if(activate_timing_adaptation){

        // With timing adaptation        

        if (Timing_Manager(0,1)==singlesupport){
 
            S_1 = round(Timing_Manager(0,0)/mpcTimeStep);
            D_1 = round(Timing_Manager(1,0)/mpcTimeStep);
            S_2 = round(Timing_Manager(2,0)/mpcTimeStep);
            D_2 = round(Timing_Manager(3,0)/mpcTimeStep);
            S_3 = round(Timing_Manager(4,0)/mpcTimeStep);
            if (floor(Timing_Manager(4,0)/mpcTimeStep) == 0){
            S_3 = 0;
            S_1 = ceil(Timing_Manager(0,0)/mpcTimeStep);
            }
            if (S_1+D_1+S_2+D_2+S_3>20 && floor(Timing_Manager(4,0)/mpcTimeStep) == 0) S_1 = S_1-1;
            if (S_1+D_1+S_2+D_2+S_3>20 && floor(Timing_Manager(4,0)/mpcTimeStep) > 0) S_3 = S_3-1;
            if (S_1+D_1+S_2+D_2+S_3<20 && floor(Timing_Manager(4,0)/mpcTimeStep) == 0) S_1 = S_1+1;
            if (S_1+D_1+S_2+D_2+S_3<20 && floor(Timing_Manager(4,0)/mpcTimeStep) > 0) S_3 = S_3+1;

            Ic.block(0,0,S_1,S_1) = Eigen::MatrixXd::Identity(S_1,S_1);
            Ic.block(S_1,S_1,D_1,D_1) = Eigen::MatrixXd::Zero(D_1,D_1);
            Ic.block(S_1+D_1,S_1+D_1,S_2,S_2) = Eigen::MatrixXd::Identity(S_2,S_2);
            Ic.block(S_1+D_1+S_2,S_1+D_1+S_2,D_2,D_2) = Eigen::MatrixXd::Zero(D_2,D_2);

            if (S_3>0) Ic.block(S_1+D_1+S_2+D_2,S_1+D_1+S_2+D_2,S_3,S_3) = Eigen::MatrixXd::Identity(S_3,S_3);

        }else{
   
            D_1 = round(Timing_Manager(0,0)/mpcTimeStep);
            S_1 = round(Timing_Manager(1,0)/mpcTimeStep);
            D_2 = round(Timing_Manager(2,0)/mpcTimeStep);
            S_2 = round(Timing_Manager(3,0)/mpcTimeStep);
            D_3 = round(Timing_Manager(4,0)/mpcTimeStep);
            if (floor(Timing_Manager(4,0)/mpcTimeStep) == 0){
            D_3 = 0;
            D_1 = ceil(Timing_Manager(0,0)/mpcTimeStep);
            }
            if (S_1+D_1+S_2+D_2+D_3>20 && floor(Timing_Manager(4,0)/mpcTimeStep) == 0) D_1 = D_1-1;
            if (S_1+D_1+S_2+D_2+D_3>20 && floor(Timing_Manager(4,0)/mpcTimeStep) > 0) D_3 = D_3-1;
            if (S_1+D_1+S_2+D_2+D_3<20 && floor(Timing_Manager(4,0)/mpcTimeStep) == 0) D_1 = D_1+1;
            if (S_1+D_1+S_2+D_2+D_3<20 && floor(Timing_Manager(4,0)/mpcTimeStep) > 0) D_3 = D_3+1;

            Ic.block(0,0,D_1,D_1) = Eigen::MatrixXd::Zero(D_1,D_1);
            Ic.block(D_1,D_1,S_1,S_1) = Eigen::MatrixXd::Identity(S_1,S_1);
            Ic.block(D_1+S_1,D_1+S_1,D_2,D_2) = Eigen::MatrixXd::Zero(D_2,D_2);
            Ic.block(D_1+S_1+D_2,D_1+S_1+D_2,S_2,S_2) = Eigen::MatrixXd::Identity(S_2,S_2);

            if (D_3>0) Ic.block(D_1+S_1+D_2+S_2,D_1+S_1+D_2+S_2,D_3,D_3) = Eigen::MatrixXd::Zero(D_3,D_3);



        }
	if((int)simulationTime/mpcTimeStep<S+D){
		Ic.block(0,0,S+D-mpcIter,S+D-mpcIter) = Eigen::MatrixXd::Zero(S+D-mpcIter,S+D-mpcIter);
	}


	Eigen::MatrixXd Cc_ = Eigen::MatrixXd::Zero(N,M);


        if (Timing_Manager(0,1) == singlesupport){

        if(Timing_Manager(5,0)==1){
        std::cout<<"A"<<std::endl;
        Cc_.block(0,0,N/2,1) = Eigen::VectorXd::Ones(N/2);
        Cc_.block(N/2,1,N/2,1) = Eigen::VectorXd::Ones(N/2);
        }else{

        Cc_.block(S_1+D_1,0,S_2+D_2,1) = Eigen::VectorXd::Ones(S_2+D_2);
        if(N-(S_1+D_1+S_2+D_2)>0) Cc_.block(S_1+D_1+S_2+D_2,1,N-(S_1+D_1+S_2+D_2),1) = Eigen::VectorXd::Ones(N-(S_1+D_1+S_2+D_2));

        }

        }else{
        Cc_.block(D_1,0,S_2+D_2,1) = Eigen::VectorXd::Ones(S_2+D_2);
        Cc_.block(D_1+S_2+D_2,1,N-(D_1+S_2+D_2),1) = Eigen::VectorXd::Ones(N-(D_1+S_2+D_2));
        }

        
        Cc = Cc_;



	rCosZmp.block(0,0,S+D-mpcIter,S+D-mpcIter) = Eigen::MatrixXd::Identity(S+D-mpcIter,S+D-mpcIter);
	rSinZmp.block(0,0,S+D-mpcIter,S+D-mpcIter) = Eigen::MatrixXd::Zero(S+D-mpcIter,S+D-mpcIter);

	for(int i=0; i<M-1; ++i){
		rCosZmp.block(S+D-mpcIter+(i*(S+D)),S+D-mpcIter+(i*(S+D)),S+D,S+D) = Eigen::MatrixXd::Identity(S+D,S+D)*cos(predictedOrientations(i+1));
		rSinZmp.block(S+D-mpcIter+(i*(S+D)),S+D-mpcIter+(i*(S+D)),S+D,S+D) = Eigen::MatrixXd::Identity(S+D,S+D)*sin(predictedOrientations(i+1));
	}

	rCosZmp.block(S+D-mpcIter+((M-1)*(S+D)),S+D-mpcIter+((M-1)*(S+D)),mpcIter,mpcIter) = Eigen::MatrixXd::Identity(S+D,S+D).block(0,0,mpcIter,mpcIter)*cos(predictedOrientations(M));
	rSinZmp.block(S+D-mpcIter+((M-1)*(S+D)),S+D-mpcIter+((M-1)*(S+D)),mpcIter,mpcIter) = Eigen::MatrixXd::Identity(S+D,S+D).block(0,0,mpcIter,mpcIter)*sin(predictedOrientations(M));


	zmpRotationMatrix << rCosZmp,rSinZmp,
						-rSinZmp,rCosZmp;





}else{  
        // No timing adaptation


	if((int)simulationTime/mpcTimeStep<S+D){
		Ic.block(0,0,S+D-mpcIter,S+D-mpcIter) = Eigen::MatrixXd::Zero(S+D-mpcIter,S+D-mpcIter);
	}
	else{
		Ic.block(0,0,S+D-mpcIter,S+D-mpcIter) = Icf.block(mpcIter,mpcIter,S+D-mpcIter,S+D-mpcIter);
	}

	for(int i=0; i<M-1; ++i){
		Ic.block(S+D-mpcIter+(i*(S+D)),S+D-mpcIter+(i*(S+D)),S+D,S+D) = Icf;
	}

	Ic.block(S+D-mpcIter+((M-1)*(S+D)),S+D-mpcIter+((M-1)*(S+D)),mpcIter,mpcIter) = Icf.block(0,0,mpcIter,mpcIter);



	for(int i=0; i<M-1; ++i){
		Cc.block(S+D-mpcIter+(i*(S+D)),i,S+D,1) = Ccf;
	}

	Cc.block(S+D-mpcIter+((M-1)*(S+D)),M-1,mpcIter,1) = Ccf.block(0,0,mpcIter,1);


	rCosZmp.block(0,0,S+D-mpcIter,S+D-mpcIter) = Eigen::MatrixXd::Identity(S+D-mpcIter,S+D-mpcIter);
	rSinZmp.block(0,0,S+D-mpcIter,S+D-mpcIter) = Eigen::MatrixXd::Zero(S+D-mpcIter,S+D-mpcIter);

	for(int i=0; i<M-1; ++i){
		rCosZmp.block(S+D-mpcIter+(i*(S+D)),S+D-mpcIter+(i*(S+D)),S+D,S+D) = Eigen::MatrixXd::Identity(S+D,S+D)*cos(predictedOrientations(i+1));
		rSinZmp.block(S+D-mpcIter+(i*(S+D)),S+D-mpcIter+(i*(S+D)),S+D,S+D) = Eigen::MatrixXd::Identity(S+D,S+D)*sin(predictedOrientations(i+1));
	}

	rCosZmp.block(S+D-mpcIter+((M-1)*(S+D)),S+D-mpcIter+((M-1)*(S+D)),mpcIter,mpcIter) = Eigen::MatrixXd::Identity(S+D,S+D).block(0,0,mpcIter,mpcIter)*cos(predictedOrientations(M));
	rSinZmp.block(S+D-mpcIter+((M-1)*(S+D)),S+D-mpcIter+((M-1)*(S+D)),mpcIter,mpcIter) = Eigen::MatrixXd::Identity(S+D,S+D).block(0,0,mpcIter,mpcIter)*sin(predictedOrientations(M));


	zmpRotationMatrix << rCosZmp,rSinZmp,
						-rSinZmp,rCosZmp;

 }       
        



        // Build constraint matrices for the QP
	AZmp.block(0,0,N,N) = Ic*P;
	AZmp.block(0,N,N,M) = -Ic*Cc;
	AZmp.block(N,N+M,N,N) = Ic*P;
	AZmp.block(N,2*N+M,N,M) = -Ic*Cc;

	AZmp = zmpRotationMatrix * AZmp;

	Eigen::VectorXd bZmpLeftTerm = Eigen::VectorXd::Zero(2*N);
	Eigen::VectorXd bZmpRightTerm = Eigen::VectorXd::Zero(2*N);

	bZmpLeftTerm << Ic*p*(footConstraintSquareWidth/2), Ic*p*(footConstraintSquareWidth/2);

	bZmpRightTerm << Ic*p*zmpPos(0), Ic*p*zmpPos(1);
	bZmpRightTerm = zmpRotationMatrix * bZmpRightTerm;

	bZmpMax =   bZmpLeftTerm - bZmpRightTerm;
	bZmpMin = - bZmpLeftTerm - bZmpRightTerm;


}

void MPCSolver::genFeasibilityConstraint(){
	AFootsteps.setZero();

	// Difference matrix (x_j - x_{j-1})
	Eigen::MatrixXd differenceMatrix = Eigen::MatrixXd::Identity(M,M);
	for (int i=0; i<M-1; ++i){
		differenceMatrix(i+1,i) = -1;
	}

	Eigen::VectorXd pFr = Eigen::VectorXd::Zero(M);
	Eigen::VectorXd pFl = Eigen::VectorXd::Zero(M);
	Eigen::VectorXd pF = Eigen::VectorXd::Ones(M);
	Eigen::MatrixXd footstepsRotationMatrix(2*M,2*M);
	Eigen::MatrixXd rCosFootsteps = Eigen::MatrixXd::Identity(M,M);
	Eigen::MatrixXd rSinFootsteps = Eigen::MatrixXd::Zero(M,M);

	AFootsteps.block(0,N,M,M) = differenceMatrix;
	AFootsteps.block(M,2*N+M,M,M) = differenceMatrix;

	rCosFootsteps(0,0) = 1;
	rSinFootsteps(0,0) = 0;

	for(int i=1; i<M; ++i){
		rCosFootsteps(i,i) = cos(predictedOrientations(i));
		rSinFootsteps(i,i) = sin(predictedOrientations(i));
	}

	footstepsRotationMatrix <<  rCosFootsteps, rSinFootsteps,
							   -rSinFootsteps, rCosFootsteps;

	AFootsteps = footstepsRotationMatrix * AFootsteps;

	if(supportFoot==true){
		for(int i=0;i<M;++i){
			if(i%2==0) pFr(i) = 1;
			else pFl(i) = 1;
		}
	}
	else{
		for(int i=0;i<M;++i){
			if(i%2==0) pFl(i) = 1;
			else pFr(i) = 1;
		}
	}

	bFootstepsMax <<  pF*deltaXMax, -pFl*deltaYIn  + pFr*deltaYOut;
	bFootstepsMin << -pF*deltaXMax, -pFl*deltaYOut + pFr*deltaYIn;
}

Eigen::VectorXd MPCSolver::solveQP() {
	int nVariables = costFunctionH.rows();
	int nConstraints = AConstraint.rows();

	qpOASES::real_t H[nVariables*nVariables];
	qpOASES::real_t g[nVariables];

	qpOASES::real_t A[nConstraints*nVariables];
	qpOASES::real_t lb[nConstraints];
	qpOASES::real_t ub[nConstraints];

	for(int i=0;i<nVariables;++i){
		for(int j=0;j<nVariables;++j){
			H[i*nVariables+j] = costFunctionH(i,j);
		}
		g[i] = costFunctionF(i);
	}

	for(int i=0;i<nConstraints;++i){
		for(int j=0;j<nVariables;++j){
			A[i*nVariables+j] = AConstraint(i,j);
		}
		lb[i] = bConstraintMin(i);
		ub[i] = bConstraintMax(i);
	}

	qpOASES::real_t xOpt[nVariables];

	qpOASES::Options options;
	options.setToMPC();
	options.printLevel=qpOASES::PL_NONE;
	qpOASES::int_t nWSR = 300;

	qp = qpOASES::QProblem(nVariables, nConstraints);
	qp.setOptions(options);
	qp.init(H,g,A,0,0,lb,ub,nWSR,NULL,NULL,NULL,NULL,NULL,NULL);

	qp.getPrimalSolution(xOpt);

	Eigen::VectorXd decisionVariables(2*(N+M));

	for(int i=0;i<2*(N+M);++i){
		decisionVariables(i) = xOpt[i];
	}

	return decisionVariables;

}

Eigen::Vector3d MPCSolver::updateState(double zmpDot, int dim, double timeStep) {
	// Update the state along the dim-th direction (0,1,2) = (x,y,z)

    double ch = cosh(omega*timeStep);
    double sh = sinh(omega*timeStep);

	Eigen::Matrix3d A_upd = Eigen::MatrixXd::Zero(3,3);
	Eigen::Vector3d B_upd = Eigen::VectorXd::Zero(3);
	A_upd<<ch,sh/omega,1-ch,omega*sh,ch,-omega*sh,0,0,1;
	B_upd<<timeStep-sh/omega,1-ch,timeStep;

    Eigen::Vector3d currentState =
    	Eigen::Vector3d( comPos(dim),comVel(dim),zmpPos(dim));

    return A_upd*currentState + B_upd*zmpDot;
}

bool MPCSolver::supportFootHasChanged(){
    if (controlIter==0 && footstepCounter>0) return true;
    else return false;
}

Eigen::VectorXd MPCSolver::getOptimalCoMPosition(){
    return comPos;
}

Eigen::VectorXd MPCSolver::getOptimalCoMVelocity(){
    return comVel;
}

Eigen::VectorXd MPCSolver::getOptimalFootsteps(){
    return predictedFootstep;
}

Eigen::VectorXd MPCSolver::getOptimalZMPPosition(){
    return zmpPos;
}

Eigen::MatrixXd MPCSolver::getPredictedZmp(){
	return predictedZmp;
}

void MPCSolver::setComTargetHeight(double height) {
	comTargetHeight = height;
}

void MPCSolver::setReferenceVelocityX(double ref) {
	vRefX = ref;
}

void MPCSolver::setReferenceVelocityY(double ref) {
	vRefY = ref;
}

void MPCSolver::setReferenceVelocityOmega(double ref) {
	omegaRef = ref;
}

void MPCSolver::setVelGain(double value) {
	qVx = value;
	qVy = value;
}

void MPCSolver::setZmpGain(double value) {
	qZ = value;
}

double MPCSolver::getOmega() {
	return omega;
}
