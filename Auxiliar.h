#pragma once
#include "Malha.h"
#include <Eigen/Dense>
//#include "random.hpp"

using Eigen::Matrix3d;
using Eigen::VectorXd;

class Auxiliar
{
public:

	long busca_linear(Malha& malha, double x, double y);
	double area_triangulo(double x1, double y1, double x2, double y2, double x3, double y3);
	double interpolar_Tri_Linear(Malha& malha, double x, double y, long elem_index, VectorXd& prop);
	bool checar_ponto_dentro_elem(unsigned long elem_index, double x, double y, Malha& malha);
};