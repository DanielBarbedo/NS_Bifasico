#pragma once
#include <string>
#include <chrono>
#include "Malha.h"
#include "Semi_Lagrangiano.h"
#include "Eigen\Sparse"
#include "Eigen\Dense"
#include "Eigen\SparseLU"

using Eigen::Matrix3d;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;
using Eigen::SparseVector;
using Eigen::Triplet;
using namespace chrono;

class NavierStokes2D
{
public:
	
	NavierStokes2D(string nome);
	void resolver_stokes_permanente();
	void resolver_stokes_transiente(double delta_t, unsigned long max_iter);
	void resolver_navier_stokes(double delta_t, unsigned long max_iter, double Re, string tag);
	
	//Projecao
	void resolver_navier_stokes_projecao(double delta_t, unsigned long max_iter, double Re, string tag);

	//Avaliação de propriedade em uma linha qualquer
	void gerar_linha_prop(int num_pontos, double xo, double yo, double x, double y);


private:

	void montar_matrizes_fixas(SparseMatrix<double>& K, SparseMatrix<double>& M, SparseMatrix<double>& Gx, SparseMatrix<double>& Gy);
	vector<Triplet<double>> montar_matriz_A(bool transiente, double delta_t, double Re);
	void aplicar_cc_A(SparseMatrix<double, Eigen::RowMajor>& A);
	void aplicar_cc_B(VectorXd& B);
	void gerar_arquivo_saida(VectorXd& B, unsigned long iter, double delta_t, double Re, string tag);

	//Avaliação de propriedade em uma linha qualquer
	void avaliar_prop(VectorXd& prop, Semi_Lagrangiano& sl, string tag, int num_pontos, double xo, double yo, double x, double y);
	
	
	//CC linha e coluna
	void aplicar_cc_linha_e_coluna(SparseMatrix<double, Eigen::RowMajor>& A);

	//Projecao
	vector<Triplet<double>> montar_matriz_B_proj(bool transiente, double delta_t, double ni);
	vector<Triplet<double>> montar_matriz_B_lumped_inversa(SparseMatrix<double, Eigen::RowMajor>& B);
	void montar_matriz_G_e_D_proj(SparseMatrix<double, Eigen::RowMajor>& G, SparseMatrix<double, Eigen::RowMajor>& D);
	void aplicar_cc_B_proj(SparseMatrix<double, Eigen::RowMajor>& B_proj);
	void aplicar_cc_D_e_G(SparseMatrix<double, Eigen::RowMajor>& D, SparseMatrix<double, Eigen::RowMajor>& G);
	void aplicar_cc_M_v(VectorXd& M_v);
	void aplicar_cc_D_v_til(VectorXd& D_v_til);
	void aplicar_cc_E(SparseMatrix<double, Eigen::RowMajor>& E);
	void aplicar_cc_D_B_inv_G(SparseMatrix<double, Eigen::RowMajor>& D_B_inv_G);
	
	Malha malha;

	unsigned long gl;
	unsigned long vert;

	SparseMatrix<double> K;
	SparseMatrix<double> M;
	SparseMatrix<double> Gx;
	SparseMatrix<double> Gy;

	//---Avaliação de propriedade---
	double linha_prop_x, linha_prop_xo, linha_prop_y, linha_prop_yo;
	int linha_prop_num_pontos;
};