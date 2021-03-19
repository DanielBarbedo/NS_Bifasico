#pragma once
#include "Malha.h"
#include "Malha_Bolha.h"
#include "NavierStokes2D.h"

class NS_Bifasico
{
public:

	NS_Bifasico(string nome, string nome_bolha);
	void resolver_navier_stokes(double delta_t, unsigned long max_iter, double mi_in, double rho_in, double mi_out, double rho_out, string tag);


private:

	//friend void NavierStokes2D::montar_matrizes_fixas(SparseMatrix<double>& K, SparseMatrix<double>& M, SparseMatrix<double>& Gx, SparseMatrix<double>& Gy);
	void calc_ts(SparseMatrix<double>& Gx, SparseMatrix<double>& Gy, VectorXd& h, VectorXd& heavyside_vec, VectorXd& nx_vec, VectorXd& ny_vec, VectorXd& tsx, VectorXd& tsy, VectorXd& kappa_x_gl, VectorXd& kappa_y_gl, VectorXd& Gx_heavy, VectorXd& Gy_heavy);
	void calc_heavyside_rho_mi_n(VectorXd& heavyside_vec, VectorXd& rho_vec, VectorXd& mi_vec, VectorXd& nx_vec, VectorXd& ny_vec, VectorXd& h, double mi_in, double rho_in, double mi_out, double rho_out);
	void montar_matrizes(SparseMatrix<double>& K, SparseMatrix<double>& M, SparseMatrix<double>& Gx, SparseMatrix<double>& Gy,
		VectorXd& heavyside_vec, VectorXd& rho_vec, VectorXd& mi_vec);
	void montar_matriz_A(SparseMatrix<double>& K, SparseMatrix<double>& M, SparseMatrix<double>& Gx, SparseMatrix<double>& Gy,
		SparseMatrix<double, Eigen::RowMajor>& A, double delta_t, double ni, double rho);
	void aplicar_cc_A(SparseMatrix<double, Eigen::RowMajor>& A);
	void aplicar_cc_B(VectorXd& B);	
	void gerar_arquivo_saida(const VectorXd& sol, const VectorXd& rho_vec, const VectorXd& mi_vec, const VectorXd& heavyside_vec, const VectorXd& tsx, const VectorXd& tsy, const VectorXd& kappa_x_gl, const VectorXd& kappa_y_gl, const VectorXd& Gx_heavy, const VectorXd& Gy_heavy, unsigned long iter, double delta_t, double ni, double rho, string tag);
	void gerar_saida_bolha(unsigned long iter, string tag);
	
	Malha malha;
	Malha_Bolha malha_bolha;

	unsigned long gl;
	unsigned long vert;
};