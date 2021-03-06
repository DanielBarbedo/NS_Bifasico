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
	double calc_heavyside(Semi_Lagrangiano& sl, unsigned long no_ref_index);
	void montar_matrizes(SparseMatrix<double>& K, SparseMatrix<double>& M, SparseMatrix<double>& Gx, SparseMatrix<double>& Gy,
		double ni1, double rho1, double ni2, double rho2, Semi_Lagrangiano& sl, VectorXd& rho_vec);
	void montar_matriz_A(SparseMatrix<double>& K, SparseMatrix<double>& M, SparseMatrix<double>& Gx, SparseMatrix<double>& Gy,
		SparseMatrix<double, Eigen::RowMajor>& A, double delta_t, double ni, double rho);
	void aplicar_cc_A(SparseMatrix<double, Eigen::RowMajor>& A);
	void aplicar_cc_B(VectorXd& B);	
	void gerar_arquivo_saida(const VectorXd& sol, const VectorXd& rho_vec, unsigned long iter, double delta_t, double ni, double rho, string tag);
	void gerar_saida_bolha(unsigned long iter, string tag);
	
	Malha malha;
	Malha_Bolha malha_bolha;

	unsigned long gl;
	unsigned long vert;
};