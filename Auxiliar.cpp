#include "Auxiliar.h"

long Auxiliar::busca_linear(Malha& malha, double x, double y)
{
	//Busca linear, percorrendo todos elementos
	for (unsigned int i = 0; i < malha.r_num_elem(); i++)
	{
		if (checar_ponto_dentro_elem(i, x, y, malha) == true) return i;
	}
	return -1;
}

double Auxiliar::area_triangulo(double x1, double y1, double x2, double y2, double x3, double y3)
{
	Matrix3d mat_A;
	mat_A << 1, x1, y1, 1, x2, y2, 1, x3, y3;

	double A = mat_A.determinant() * 0.5;
	return abs(A);
}

double Auxiliar::interpolar_Tri_Linear(Malha& malha, double x, double y,
	long elem_index, VectorXd& prop)
{
	Elemento elem = malha.r_elem(elem_index);

	No no1 = malha.r_no(elem.nos[0]);
	No no2 = malha.r_no(elem.nos[1]);
	No no3 = malha.r_no(elem.nos[2]);

	double peso1 = ((no2.y - no3.y) * (x - no3.x) + (no3.x - no2.x) * (y - no3.y)) /
		((no2.y - no3.y) * (no1.x - no3.x) + (no3.x - no2.x) * (no1.y - no3.y));
	double peso2 = ((no3.y - no1.y) * (x - no3.x) + (no1.x - no3.x) * (y - no3.y)) /
		((no2.y - no3.y) * (no1.x - no3.x) + (no3.x - no2.x) * (no1.y - no3.y));
	double peso3 = 1 - peso1 - peso2;

	return prop[no1.id] * peso1 + prop[no2.id] * peso2 + prop[no3.id] * peso3;
}

bool Auxiliar::checar_ponto_dentro_elem(unsigned long elem_index, double x, double y, Malha& malha)
{
	Elemento elem = malha.r_elem(elem_index);

	No no1 = malha.r_no(elem.nos[0]);
	No no2 = malha.r_no(elem.nos[1]);
	No no3 = malha.r_no(elem.nos[2]);

	double A = area_triangulo(no1.x, no1.y, no2.x, no2.y, no3.x, no3.y);
	double A1 = area_triangulo(x, y, no2.x, no2.y, no3.x, no3.y);
	double A2 = area_triangulo(no1.x, no1.y, x, y, no3.x, no3.y);
	double A3 = area_triangulo(no1.x, no1.y, no2.x, no2.y, x, y);

	double soma = A1 + A2 + A3;
	double percentual = soma / A;

	if (percentual >= 0.999999 && percentual <= 1.00001) return true;
	else return false;
}