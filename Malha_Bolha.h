#pragma once
#include <string>
#include <fstream>
#include <iostream>
#include "Malha.h"
#include "Eigen\Dense"

using namespace std;
using namespace Eigen;

class Malha_Bolha
{
public:
	void ler_malha(string nome_arquivo);
	void mostrar_malha();
	void atualizar_posicao(const VectorXd& vx, const VectorXd& vy, double dt);
	unsigned long r_num_nos() { return (unsigned long)m_no_vec.size(); }
	No r_no(unsigned long no_index) { return m_no_vec[no_index]; }
	No r_no_centro() { return m_no_centro; }
	Elemento r_elem(unsigned long elem_index) { return m_elem_vec[elem_index]; }

private:

	void ler_nos(ifstream& arquivo);
	void ler_elementos(ifstream& arquivo);	

	No m_no_centro;
	vector<No> m_no_vec;
	vector<Elemento> m_elem_vec;
};