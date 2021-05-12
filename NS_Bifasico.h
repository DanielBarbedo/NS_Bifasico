#pragma once
#include "Malha.h"
#include "Malha_Bolha.h"
#include "NavierStokes2D.h"

class NS_Bifasico
{
public:

	NS_Bifasico(string nome, string nome_bolha);
	void resolver_navier_stokes(double delta_t, unsigned long max_iter, double Reynolds, double mi_in, double rho_in, double mi_out, double rho_out, string tag);

	//Avaliação de propriedade em uma linha qualquer
	void gerar_linha_prop(int num_pontos, double xo, double yo, double x, double y);

private:

	//friend void NavierStokes2D::montar_matrizes_fixas(SparseMatrix<double>& K, SparseMatrix<double>& M, SparseMatrix<double>& Gx, SparseMatrix<double>& Gy);
	void calc_ts(SparseMatrix<double>& Gx, SparseMatrix<double>& Gy, VectorXd& h, VectorXd& heavyside_vec, VectorXd& nx_vec, VectorXd& ny_vec, VectorXd& tsx, VectorXd& tsy, VectorXd& kappa_gl, VectorXd& Gx_heavy, VectorXd& Gy_heavy);
	void calc_heavyside_rho_mi_n(VectorXd& heavyside_vec, VectorXd& rho_vec, VectorXd& mi_vec, VectorXd& nx_vec, VectorXd& ny_vec, VectorXd& h, double mi_in, double rho_in, double mi_out, double rho_out);
	void montar_matrizes(SparseMatrix<double>& Kx, SparseMatrix<double>& Ky, SparseMatrix<double>& Kxy, SparseMatrix<double>& M, SparseMatrix<double>& Gx, SparseMatrix<double>& Gy, VectorXd& rho_vec, VectorXd& mi_vec);
	void montar_matriz_A(SparseMatrix<double>& Kx, SparseMatrix<double>& Ky, SparseMatrix<double>& Kxy, SparseMatrix<double>& M, SparseMatrix<double>& Gx, SparseMatrix<double>& Gy, SparseMatrix<double, Eigen::RowMajor>& A, double delta_t, double Reynolds);
	void aplicar_cc_A(SparseMatrix<double, Eigen::RowMajor>& A);
	void aplicar_cc_B(VectorXd& B);	
	void gerar_arquivo_saida(const VectorXd& sol, const VectorXd& rho_vec, const VectorXd& mi_vec, const VectorXd& heavyside_vec, const VectorXd& tsx, const VectorXd& tsy, const VectorXd& kappa_gl, const VectorXd& Gx_heavy, const VectorXd& Gy_heavy, unsigned long iter, double delta_t, double Reynolds, string tag);
	void gerar_saida_bolha(unsigned long iter, double delta_t, double Reynolds, string tag);
	void plotar_diametros(const VectorXd& vx, const VectorXd& vy, string tag, unsigned long iter, double dt);

	//Avaliação de propriedade em uma linha qualquer
	void avaliar_prop(VectorXd& prop, Semi_Lagrangiano& sl, string tag, int num_pontos, double xo, double yo, double x, double y);
	
	Malha malha;
	Malha_Bolha malha_bolha;

	unsigned long gl;
	unsigned long vert;

	//---Avaliação de propriedade---
	double linha_prop_x, linha_prop_xo, linha_prop_y, linha_prop_yo;
	int linha_prop_num_pontos;
};