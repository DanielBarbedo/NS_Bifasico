#include "Malha.h"

void Malha::ler_malha(string nome_arquivo)
{
	m_nome = nome_arquivo;
	m_num_vertices = 0;

	m_tags_cc.inflow = -1;
	m_tags_cc.outflow = -1;
	m_tags_cc.wall = -1;
	m_tags_cc.inflow_vert = -1;
	m_tags_cc.outflow_vert = -1;
	m_tags_cc.wall_slip_ver = -1;
	m_tags_cc.wall_slip_hor = -1;

	ifstream arquivo;
	arquivo.open(nome_arquivo + ".msh");	
	if (arquivo.is_open() == false) std::cout << "Abertura de arquivo da malha falhou" << endl;

	string str;

	while (true)
	{
		arquivo >> str;
		if (str == "$PhysicalNames") ler_tags_cond_contorno(arquivo);
		if (str == "$Nodes") ler_nos(arquivo);
		if (str == "$Elements") ler_elementos(arquivo);
		if (arquivo.eof() == true) break;
	}
	arquivo.close();
}

void Malha::ler_tags_cond_contorno(ifstream& arquivo) //modificada
{
	//Obtem-se e descarta-se o número de tags. Só são importantes
	//as tags que descrevem uma condicao de contorno
	string str;
	arquivo >> str;

	int tag;
	//A partir deste ponto, as linhas são estruturadas em 3 elementos:
	//A classificação da "entidade", sua tag, e a descrição
	while (true)
	{
		//Pega-se o primeiro numero e descarta-se, este numero representa
		//o tipo do grupo, se são pontos, segmentos de reta, superfícies
		//ou volumes
		arquivo >> str;
		if (str == "$EndPhysicalNames") break;

		//O segundo número é a tag para esta categoria.
		arquivo >> str;
		tag = stoul(str);

		//A terceira palavra é a descrição. Só é gravada a tag referente
		//ao inflow, outflow ou wall
		arquivo >> str;

		cout << str;

		if (str == "\"inflow\"") m_tags_cc.inflow = tag;
		else if (str == "\"outflow\"") m_tags_cc.outflow = tag;
		else if (str == "\"inflow_vert\"") m_tags_cc.inflow_vert = tag;
		else if (str == "\"outflow_vert\"") m_tags_cc.outflow_vert = tag;
		else if (str == "\"wall\"") m_tags_cc.wall = tag;
		else if (str == "\"wall_slip_ver\"") m_tags_cc.wall_slip_ver = tag;
		else if (str == "\"wall_slip_hor\"") m_tags_cc.wall_slip_hor = tag;
		else if (str == "\"corner_lid\"") m_tags_cc.corner_lid = tag;
	}

	//if (m_tags_cc.inflow == -1 || m_tags_cc.outflow == -1 || m_tags_cc.wall == -1)
	//	cout << "Erro na leitura das tags das condicoes de contorno";

	cout << "Mostrando as tags e valores obtidos" << endl;
	cout << "Inflow: " << m_tags_cc.inflow << endl;
	cout << "Outflow: " << m_tags_cc.outflow << endl;
	cout << "Inflow for a vertical flow: " << m_tags_cc.inflow_vert << endl;
	cout << "Outflow for a vertical flow: " << m_tags_cc.outflow_vert << endl;
	cout << "Wall: " << m_tags_cc.wall << endl;
	cout << "Horizontal Wall Slip: " << m_tags_cc.wall_slip_hor << endl;
	cout << "Vertical Wall Slip: " << m_tags_cc.wall_slip_ver << endl;
	cout << "Lid driven null pressure corner: " << m_tags_cc.corner_lid << endl;
}

void Malha::ler_nos(ifstream& arquivo)
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

void Malha::ler_elementos(ifstream& arquivo)
{
	string str;
	unordered_map<unsigned long, int> cc_geral_map;
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

		if (tipo == 15)
		{
			//Descarta-se o proximo número, que é a quantidade de tags (2)
			arquivo >> str;			
			//Este número é a tag do elemento, que será comparada com as tags
			//das condicoes de contorno obtidas em outra função e armazenadas
			arquivo >> str;
			int tag = stoi(str);
			//Proximo numero é descartado
			arquivo >> str;
			arquivo >> str;
			unsigned long no_index = stoul(str);
			no_index = no_index - 1; //Correção da numeração. Queremos começar no nó 0 ao invés do 1.
			checar_cond_contorno(this->r_no(no_index), tag, cc_geral_map);
		}

		//1 -> elemento de linha 2 nós
		//8 -> elemento de linha 3 nós
		//Os elementos de linha carregam consigo tags das condicoes de contorno
		if (tipo == 1 || tipo == 8)
		{
			//Descarta-se o proximo número, que é a quantidade de tags (2)
			arquivo >> str;

			//Este número é a tag do elemento, que será comparada com as tags
			//das condicoes de contorno obtidas em outra função e armazenadas
			arquivo >> str;
			int tag = stoi(str);
			
			//Proximo numero é descartado
			arquivo >> str;

			//Os próximos numeros são os nós do elemento de linha
			int num_nos = 0;
			if (tipo == 1) num_nos = 2;
			if (tipo == 8) num_nos = 3;
			
			for (int i = 0; i < num_nos; i++)
			{
				arquivo >> str;
				unsigned long no_index = stoul(str);
				no_index = no_index - 1; //Correção da numeração. Queremos começar no nó 0 ao invés do 1.
				checar_cond_contorno(this->r_no(no_index), tag, cc_geral_map);

				//Armazena-se os nós de todos os elementos de linha. Estes nós definem o contorno
				//completo da malha, e são necessários para aplicação de algumas cond. contorno
				bool no_duplicado = false;
				for (unsigned int i = 0; i < m_no_contorno_vec.size(); i++)
				{
					if (m_no_contorno_vec[i] == no_index) no_duplicado = true;
				}
				if (no_duplicado == false) m_no_contorno_vec.push_back(no_index);
			}
		}
		
		//2 -> elemento triangular linear, 3 nós
		//9 -> elemento triangular quadrático, 6 nós
		if (tipo == 2 || tipo == 9)
		{
			struct Elemento elem;
			int num_elem = 0;
			if (tipo == 2)
			{
				num_elem = 3;
				elem.tipo = Tipo_Elem::Tri_Linear;
			}
			if (tipo == 9)
			{
				num_elem = 6;
				elem.tipo = Tipo_Elem::Tri_Quadr;
			}
			
			for (int i = 0; i < 3; i++) arquivo >> str;
			for (int i = 0; i < num_elem; i++)
			{
				arquivo >> str;
				elem.nos.push_back(stoul(str) - 1);
			}
			m_elem_vec.push_back(elem);
		}
	}

	//Finalizando a leitura dos elementos, faz-se a separação das 
	//condições de contorno, de acordo com a propriedade apropriada.
	for (unordered_map<unsigned long, int>::iterator it = cc_geral_map.begin(); it != cc_geral_map.end(); it++)
	{
		if (it->second == m_tags_cc.inflow)
		{
			m_cc_vx_map.insert(pair<unsigned long, double>(it->first, 1));
			m_cc_vy_map.insert(pair<unsigned long, double>(it->first, 0));
		}
		else if (it->second == m_tags_cc.inflow_vert)
		{
			m_cc_vx_map.insert(pair<unsigned long, double>(it->first, 0));
			m_cc_vy_map.insert(pair<unsigned long, double>(it->first, 1));
		}
		else if (it->second == m_tags_cc.outflow)
		{
			m_cc_p_map.insert(pair<unsigned long, double>(it->first, 0));
		}
		else if (it->second == m_tags_cc.outflow_vert)
		{
			m_cc_p_map.insert(pair<unsigned long, double>(it->first, 0));
		}
		else if (it->second == m_tags_cc.wall)
		{
			m_cc_vx_map.insert(pair<unsigned long, double>(it->first, 0));
			m_cc_vy_map.insert(pair<unsigned long, double>(it->first, 0));
		}
		else if (it->second == m_tags_cc.wall_slip_hor)
		{
			m_cc_vx_map.insert(pair<unsigned long, double>(it->first, 1));
			m_cc_vy_map.insert(pair<unsigned long, double>(it->first, 0));
		}
		else if (it->second == m_tags_cc.wall_slip_ver)
		{
			m_cc_vx_map.insert(pair<unsigned long, double>(it->first, 0));
			m_cc_vy_map.insert(pair<unsigned long, double>(it->first, 1));
		}
		else if (it->second == m_tags_cc.corner_lid)
		{
			m_cc_vx_map.insert(pair<unsigned long, double>(it->first, 0));
			m_cc_vy_map.insert(pair<unsigned long, double>(it->first, 0));
			m_cc_p_map.insert(pair<unsigned long, double>(it->first, 0));
		}
		else cout << "Erro na construcao das ccs por variavel" << endl;
	}
}

void Malha::checar_cond_contorno(No no, int tag, unordered_map<unsigned long, int>& cc_geral_map)
{
	//Tenta-se inserir o elemento no mapa geral		
	pair<unordered_map<unsigned long, int>::iterator, bool> par = cc_geral_map.insert(pair<unsigned long, int>(no.id, tag));
	//Ele retorna um par<iterador, bool>. Se o bool for true, significa que o elemento
	//era único e foi inserido. Retornamos.
	if (par.second == true) return;
	//Caso contrário o elemento já existe. Neste caso tem-se dois cenários
	//1 - O elemento a ser adicionado é um elemento inflow ou outflow. Neste caso
	//é um elemento duplicado (adicionado antes)
	//2 - O elemento a ser adicionado é um elemento wall. Neste caso, ele tanto
	//pode ser uma duplicata, como pode ser um caso de intersecção de condição
	//de contorno. No caso da interseção com inflow, queremos que wall tome prioridade
	else if (tag == m_tags_cc.wall &&
		(par.first->second == m_tags_cc.inflow || par.first->second == m_tags_cc.outflow
			|| par.first->second == m_tags_cc.outflow_vert || par.first->second == m_tags_cc.inflow_vert))
		par.first->second = m_tags_cc.wall;
	else if (tag == m_tags_cc.corner_lid) par.first->second = m_tags_cc.corner_lid;
	//else if (tag == m_tags_cc.wall_slip_hor) par.first->second = m_tags_cc.wall_slip_hor;
	//else if (tag == m_tags_cc.wall_slip_ver) par.first->second = m_tags_cc.wall_slip_ver;
}	

void Malha::mostrar_malha()
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

	cout << "Existem " << m_cc_vx_map.size() << " nos com condicoes de contorno vx definidas" << endl;
	for (unordered_map<unsigned long, double>::iterator it = m_cc_vx_map.begin(); it != m_cc_vx_map.end(); it++)
	{
		cout << "No e valor da condicao de contorno: " << it->first << " " << it->second << " " << endl;
	}
	cout << endl;

	cout << "Existem " << m_cc_vy_map.size() << " nos com condicoes de contorno vy definidas" << endl;
	for (unordered_map<unsigned long, double>::iterator it = m_cc_vy_map.begin(); it != m_cc_vy_map.end(); it++)
	{
		cout << "No e valor da condicao de contorno: " << it->first << " " << it->second << " " << endl;
	}
	cout << endl;

	cout << "Existem " << m_cc_p_map.size() << " nos com condicoes de contorno p definidas" << endl;
	for (unordered_map<unsigned long, double>::iterator it = m_cc_p_map.begin(); it != m_cc_p_map.end(); it++)
	{
		cout << "No e valor da condicao de contorno: " << it->first << " " << it->second << " " << endl;
	}
	cout << endl;

	cout << "Existem " << m_no_contorno_vec.size() << " nos no contorno da malha" << endl;
	for (unsigned long i = 0; i < m_no_contorno_vec.size(); i++)
	{
		cout << "No no contorno da malha: " << m_no_contorno_vec[i] << endl;
	}
	cout << endl;
}

Elemento Malha::r_elem(int index)
{
	if (index == -1) cout << "Malha::r_elem, indice invalido" << endl;
	return m_elem_vec[index];
}

No Malha::r_no(unsigned long no_ID)
{
	return m_no_vec[no_ID];
}

void Malha::gerar_elemento_mini()
{
	m_num_vertices = m_no_vec.size();
	
	for(unsigned int i = 0; i < m_elem_vec.size(); i++)
	{
		No no1 = r_no(m_elem_vec[i].nos[0]);
		No no2 = r_no(m_elem_vec[i].nos[1]);
		No no3 = r_no(m_elem_vec[i].nos[2]);

		No novo;		
		novo.id = m_no_vec.size();
		novo.x = (no1.x + no2.x + no3.x) / 3;
		novo.y = (no1.y + no2.y + no3.y) / 3;
		novo.z = (no1.z + no2.z + no3.z) / 3;

		m_elem_vec[i].nos.push_back(novo.id);

		m_no_vec.push_back(novo);		
	}
}