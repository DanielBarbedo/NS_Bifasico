#pragma once
#include "Malha.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unordered_map>
#include <chrono>

using namespace std;
using Eigen::Matrix3d;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;
using Eigen::Triplet;
using namespace chrono;

struct Elem_Oposto
{
	unsigned long no;
	long elem;
};

class Semi_Lagrangiano
{
public:
	
	Semi_Lagrangiano(Malha& malha);
	void semi_lagrangiano(Malha& malha, VectorXd& vx, VectorXd& vy,	double delta_t);
	void multiplicar_SL(VectorXd& prop);
	
	//Estão aqui pq precisamos interpolar a velocidade pra fazer a malha da interface se mover
	long busca_linear(Malha& malha, double x, double y);

private:

	void montar_Auxiliar(Malha& malha);
	void mostrar_Auxiliar();
	bool checar_oposicao(Elemento elem_ref, Elemento elem, int a, int b);

	bool checar_contorno(Malha& malha, Elemento elem, int a, int b);
	void achar_novo_valor(unsigned long index_origem, Malha& malha, double x, double y);
	void interpolar_Tri_Linear(Malha& malha, double x, double y, long index_atual, long index_origem);
	void interpolar_contorno(Malha& malha, long index_atual, long index_origem, double xo, double yo, double x, double y);

	//unsigned long buscar_no_prox2(Malha& malha, double x, double y);

	//long busca_linear(Malha& malha, double x, double y);
	double area_triangulo(double x1, double y1, double x2, double y2, double x3, double y3);	
	bool checar_ponto_dentro_elem(unsigned long elem_index, double x, double y, Malha& malha);

	vector< vector<unsigned long> > nos_elem_vec;
	vector< vector<Elem_Oposto> > elem_opostos;

	SparseMatrix<double> mat_SL;
	vector<Triplet<double>> SL_triplets;

	unsigned long contagem_busca_linear;
	float duracao_busca_linear;
};