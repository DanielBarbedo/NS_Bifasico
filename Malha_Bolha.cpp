#include "Malha_Bolha.h"

void Malha_Bolha::ler_malha(string nome_arquivo)
{
	//m_nome = nome_arquivo;	
	ifstream arquivo;
	arquivo.open(nome_arquivo + ".msh");
	if (arquivo.is_open() == false) std::cout << "Abertura de arquivo da malha falhou" << endl;

	string str;

	while (true)
	{
		arquivo >> str;
		//if (str == "$PhysicalNames") ler_tags_cond_contorno(arquivo);
		if (str == "$Nodes") ler_nos(arquivo);
		if (str == "$Elements") ler_elementos(arquivo);
		if (arquivo.eof() == true) break;
	}
	arquivo.close();
}

void Malha_Bolha::ler_nos(ifstream& arquivo)
{
	string str;
	struct No no;

	arquivo >> str;
	int num_nos = stoi(str);

	for (int i = 0; i < num_nos; i++)
	{
		arquivo >> str;
		no.id = stoul(str);
		no.id = no.id - 1;
		arquivo >> str;
		no.x = stod(str);
		arquivo >> str;
		no.y = stod(str);
		arquivo >> str;
		no.z = stod(str);

		m_no_vec.push_back(no);
	}
}

void Malha_Bolha::ler_elementos(ifstream& arquivo)
{
	string str;
	arquivo >> str;	
	int num_elementos = stoi(str);

	for (int j = 0; j < num_elementos; j++)
	{
		//O primeiro passo é extrair o numero do elemento, o seu id.
		//A contagem realizada pelo gerador de malha inclui pontos e segmentos
		//de reta, então o id não é útil para o trabalho do software
		arquivo >> str;

		//O numero seguinte nos dá o tipo do elemento, de acordo com o seu tipo, descobre-se
		//quantos números restantes (tags e coordenadas) existem em uma linha
		//O número de classificações é constante, mas o número de coordenadas muda.
		arquivo >> str;
		int tipo;
		tipo = stoul(str);

		if (tipo == 15) //Ponto. Não é útil para nada
		{
			//Descarta-se o proximo número, que é a quantidade de tags (2)
			arquivo >> str;			
			//Este número é a tag do elemento, que será comparada com as tags
			//das condicoes de contorno obtidas em outra função e armazenadas
			arquivo >> str;
			//Proximo numero é descartado
			arquivo >> str;
			arquivo >> str;
			m_no_centro = r_no(stoul(str));
		}

		//1 -> elemento de linha 2 nós
		//8 -> elemento de linha 3 nós
		//Os elementos de linha carregam consigo tags das condicoes de contorno
		if (tipo == 1 || tipo == 8)
		{
			//Descarta-se o proximo número, que é a quantidade de tags (2)
			arquivo >> str;
			//Este número é a tag do elemento, não tem utilidade por hora
			arquivo >> str;			
			//Proximo numero é descartado
			arquivo >> str;

			//Os próximos numeros são os nós do elemento de linha
			int num_nos = 0;
			if (tipo == 1) num_nos = 2;
			if (tipo == 8) num_nos = 3;
			
			struct Elemento elem;
			elem.tipo = Tipo_Elem::Linha;

			for (int i = 0; i < num_nos; i++)
			{			
				arquivo >> str;
				elem.nos.push_back(stoul(str) - 1);
			}

			m_elem_vec.push_back(elem);
		}
	}
}

void Malha_Bolha::mostrar_malha()
{
	cout << "Existem " << m_no_vec.size() << " nos registrados" << endl;
	for (unsigned int i = 0; i < m_no_vec.size(); i++)
	{
		cout << "No id, x, y, z:" << m_no_vec[i].id << " ";
		cout << m_no_vec[i].x << " " << m_no_vec[i].y << " " << m_no_vec[i].z << endl;
	}
	cout << endl;

	cout << "Existem " << m_elem_vec.size() << " elementos registrados" << endl;
	for (unsigned int i = 0; i < m_elem_vec.size(); i++)
	{
		cout << "Elemento id, e seus nos: " << i << ": ";
		for (unsigned int j = 0; j < m_elem_vec[i].nos.size(); j++)	cout << m_elem_vec[i].nos[j] << " ";
		cout << endl;
	}
	cout << endl;
}

void Malha_Bolha::atualizar_posicao(Malha& malha, Semi_Lagrangiano& sl, const VectorXd& vx, const VectorXd& vy, double dt)
{	
	for (unsigned int i = 0; i < m_no_vec.size(); i++)
	{
		long el_index = sl.busca_linear(malha, m_no_vec[i].x, m_no_vec[i].y);
		if (el_index == -1) cout << "Malha_Bolha::atualizar_posicao, elemento fora do dominio" << endl;
		double vx_ponto = interpolar_Tri_Linear(malha, m_no_vec[i].x, m_no_vec[i].y, el_index, vx);
		double vy_ponto = interpolar_Tri_Linear(malha, m_no_vec[i].x, m_no_vec[i].y, el_index, vy);
		m_no_vec[i].x += vx_ponto * dt;
		m_no_vec[i].y += vy_ponto * dt;
	}
}

double Malha_Bolha::interpolar_Tri_Linear(Malha& malha, double x, double y, long elem_index, const VectorXd& prop)
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