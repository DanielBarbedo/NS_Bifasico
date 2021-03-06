#pragma once
#include <vector>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

enum class Tipo_Elem
{
	Linha,
	Tri_Linear,
	Tri_Quadr,
	Mini
};

struct No
{
	unsigned long id;
	double x;
	double y;
	double z;
};

struct Elemento
{
	enum Tipo_Elem tipo;
	vector<unsigned long> nos;
};

struct Tags_CC
{
	int inflow;
	int outflow;
	int inflow_vert;
	int outflow_vert;
	int wall;
	int wall_slip_ver;
	int wall_slip_hor;	
	int corner_lid;
};

class Malha
{
public:

	void ler_malha(string nome_arquivo);
	string r_nome() { return m_nome; }
	void mostrar_malha();
	unsigned long r_num_elem() { return (unsigned long)m_elem_vec.size(); }
	unsigned long r_num_nos() { return (unsigned long)m_no_vec.size(); }
	unsigned long r_num_vertices() { return m_num_vertices; }
	Elemento r_elem(int index);
	No r_no(unsigned long no_ID);
	unordered_map<unsigned long, double> r_cc_vx() { return m_cc_vx_map; }
	unordered_map<unsigned long, double> r_cc_vy() { return m_cc_vy_map; }
	unordered_map<unsigned long, double> r_cc_p() { return m_cc_p_map; }
	vector<unsigned long> r_no_contorno() { return m_no_contorno_vec; }

	void gerar_elemento_mini();

private:
		
	void ler_nos(ifstream& arquivo);
	void ler_elementos(ifstream& arquivo);
	void ler_tags_cond_contorno(ifstream& arquivo);
	void checar_cond_contorno(No no, int tag, unordered_map<unsigned long, int>& cc_geral_map);

	unsigned long m_num_vertices;
	
	string m_nome;
	vector<No> m_no_vec;	
	vector<Elemento> m_elem_vec;
	vector<unsigned long> m_no_contorno_vec;
	Tags_CC m_tags_cc;
	
	unordered_map<unsigned long, double> m_cc_vx_map;
	unordered_map<unsigned long, double> m_cc_vy_map;
	unordered_map<unsigned long, double> m_cc_p_map;
};