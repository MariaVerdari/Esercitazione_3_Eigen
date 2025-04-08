#include <iostream>
#include "Eigen/Eigen"
#include <iomanip>

using namespace std;
using namespace Eigen;

int main()
{
	cout<<endl;
	
	Vector2d real_x(-1.0e+0, -1.0e+00);
	double real_norm = real_x.norm();

	
	// Primo sistema
	Matrix2d arrMat[3];
	Vector2d arrVec[3];
	
	Matrix2d  A1 = Matrix2d::Zero();
	A1 << 5.547001962252291e-01, -3.770900990025203e-02, 8.320502943378437e-01,
	-9.992887623566787e-01;
	arrMat[0] = A1;

	Vector2d b1(-5.169911863249772e-01, 1.672384680188350e-01);
	arrVec[0] = b1;
	// Secondo sistema
	
	Matrix2d  A2 = Matrix2d::Zero();
	A2 << 5.547001962252291e-01, -5.540607316466765e-01, 8.320502943378437e-01,
	-8.324762492991313e-01;
	arrMat[1] = A2;
	Vector2d b2(-6.394645785530173e-04, 4.259549612877223e-04);
	
	arrVec[1] = b2;
	
		// Terzo sistema
	Matrix2d  A3 = Matrix2d::Zero();
	A3 << 5.547001962252291e-01, -5.547001955851905e-01, 8.320502943378437e-01,
	-8.320502947645361e-01;
	arrMat[2] = A3;
	Vector2d b3(-6.400391328043042e-10, 4.266924591433963e-10);
	arrVec[2] = b3;
	
	std::string arrStr[3] = {"primo", "secondo", "terzo"};
	
	Vector2d arrLU[3];
	Vector2d arrDLU[3];
	double arrELU[3];
	Vector2d arrQR[3];
	Vector2d arrDQR[3];
	double arrEQR[3];
	
	

	for(int i=0; i<3; i++) {
		arrLU[i]= arrMat[i].lu().solve(arrVec[i]); // con pivoting parziale
		cout << setprecision(16)<< scientific<< "La soluzione del " <<arrStr[i] <<" sistema lineare calcolata con la fattorizzazione PALU è: "<< endl<< arrLU[i]<<endl<<endl;
		arrDLU[i] = real_x - arrLU[i];
		arrELU[i]= arrDLU[i].norm();
		cout << "L'errore relativo del "<<arrStr[i]  <<" sistema lineare con la fattorizzazione PALU è: "<< endl<< arrELU[i]/real_norm <<endl<<endl;

		
		arrQR[i]= arrMat[i].householderQr().solve(arrVec[i]); 
		cout << "La soluzione del " <<arrStr[i] <<" sistema lineare calcolata con la fattorizzazione QR è: "<< endl<< arrQR[i]<<endl<<endl;
		arrDQR[i] = real_x - arrQR[i];
		arrEQR[i]= arrDQR[i].norm();
		cout << "L'errore relativo del "<<arrStr[i]  <<" sistema lineare con la fattorizzazione QR è: "<< endl<< arrEQR[i]/real_norm<<endl<<endl;

		
	}
	
	

	

    return 0;
}