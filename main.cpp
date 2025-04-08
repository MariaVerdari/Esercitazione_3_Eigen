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
	Matrix2d  A1 = Matrix2d::Zero();
	A1 << 5.547001962252291e-01, -3.770900990025203e-02, 8.320502943378437e-01,
	-9.992887623566787e-01;
	Vector2d b1(-5.169911863249772e-01, 1.672384680188350e-01);
	
	Vector2d x1_lu = A1.lu().solve(b1); // con pivoting parziale
	cout << setprecision(16)<< scientific<< "La soluzione del primo sistema lineare calcolata con la fattorizzazione PALU è: "<< endl<< x1_lu<<endl<<endl;
	
	Vector2d  diff1_lu = real_x -x1_lu;
	double err1_lu = diff1_lu.norm()/real_norm;
	
	cout << "L'errore relativo del primo sistema lineare con la fattorizzazione PALU è: "<< endl<< err1_lu<<endl<<endl;
	

	Vector2d x1_qr = A1.householderQr().solve(b1);
	cout << "La soluzione del primo sistema lineare calcolata con la fattorizzazione QR è: "<< endl<< x1_qr<<endl<<endl;
	
	Vector2d  diff1_qr = real_x -x1_qr;
	double err1_qr = diff1_qr.norm()/real_norm;
	
	cout << "L'errore relativo del primo sistema lineare con la fattorizzazione QR è: "<< endl<< err1_qr<<endl<<endl<<endl;
	
	// Secondo sistema
	Matrix2d  A2 = Matrix2d::Zero();
	A2 << 5.547001962252291e-01, -5.540607316466765e-01, 8.320502943378437e-01,
	-8.324762492991313e-01;
	Vector2d b2(-6.394645785530173e-04, 4.259549612877223e-04);
	
	Vector2d x2_lu = A2.lu().solve(b2); // con pivoting parziale
	cout << "La soluzione del secondo sistema lineare calcolata con la fattorizzazione PALU è: "<< endl<< x2_lu<<endl<<endl;
	
	Vector2d  diff2_lu = real_x -x2_lu;
	double err2_lu = diff2_lu.norm()/real_norm;
	
	cout << "L'errore relativo del secondo sistema lineare con la fattorizzazione PALU è: "<< endl<< err2_lu<<endl<<endl;
	
	
	
	Vector2d x2_qr = A2.householderQr().solve(b2);
	cout << "La soluzione del secondo sistema lineare calcolata con la fattorizzazione QR è: "<< endl<< x2_qr<<endl<<endl;
	
	Vector2d  diff2_qr = real_x -x2_qr;
	double err2_qr = diff2_qr.norm()/real_norm;
	
	cout << "L'errore relativo del secondo sistema lineare con la fattorizzazione QR è: "<< endl<< err2_qr<<endl<<endl<<endl;
	
	
	// Terzo sistema
	Matrix2d  A3 = Matrix2d::Zero();
	A3 << 5.547001962252291e-01, -5.547001955851905e-01, 8.320502943378437e-01,
	-8.320502947645361e-01;
	Vector2d b3(-6.400391328043042e-10, 4.266924591433963e-10);
	
	Vector2d x3_lu = A3.lu().solve(b3); // con pivoting parziale
	cout << "La soluzione del terzo sistema lineare calcolata con la fattorizzazione PALU è: "<< endl<< x3_lu<<endl<<endl;
	
	Vector2d  diff3_lu = real_x -x3_lu;
	double err3_lu = diff3_lu.norm()/real_norm;
	
	cout << "L'errore relativo del terzo sistema lineare con la fattorizzazione PALU è: "<< endl<< err3_lu<<endl<<endl;
	
	
	
	Vector2d x3_qr = A3.householderQr().solve(b3);
	cout << "La soluzione del terzo sistema lineare calcolata con la fattorizzazione QR è: "<< endl<< x3_qr<<endl<<endl;
	
	Vector2d  diff3_qr = real_x -x3_qr;
	double err3_qr = diff3_qr.norm()/real_norm;
	
	cout << "L'errore relativo del terzo sistema lineare con la fattorizzazione QR è: "<< endl<< err3_qr<<endl<<endl<<endl;

	

    return 0;
}
 
 

