#include "Semi_Lagrangiano.h"

Semi_Lagrangiano::Semi_Lagrangiano(Malha& malha)
{
	montar_Auxiliar(malha);
	//mostrar_Auxiliar();
}

void Semi_Lagrangiano::montar_Auxiliar(Malha& malha)
{
	//Aqui fazemos uma lista dos elementos que possuem o nó como parte
	//elem_nos, é uma lista, cujos índices são os índices dos nós e o conteúdo
	//é uma lista dos elementos que possuem o nó como composição

	for (unsigned long i = 0; i < malha.r_num_nos(); i++)
	{
		vector<unsigned long> elem_vec;
		nos_elem_vec.push_back(elem_vec);
	}	
	for (unsigned int j = 0; j < malha.r_num_elem(); j++)
	{
		Elemento el = malha.r_elem(j);
		nos_elem_vec[el.nos[0]].push_back(j);
		nos_elem_vec[el.nos[1]].push_back(j);
		nos_elem_vec[el.nos[2]].push_back(j);
	}	
	unsigned long num = 0;
	for (unsigned int i = malha.r_num_vertices(); i < malha.r_num_nos(); i++)
	{
		nos_elem_vec[i].push_back(num);
		num++;
	}

	//	//Aqui construimos uma lista dos elementos opostos
	//	//Tomamos um elemento como referência, checamos seus nós e os elementos
	//	//opostos a estes nós. Caso o nó não tenha um elemento oposto
	//	//(está no contorno da malha), utiliza-se o -1
	
	for (unsigned int j = 0; j < malha.r_num_elem(); j++)
	{
		Elemento el_ref = malha.r_elem(j);

		vector<unsigned long> el_no_comp;
		el_no_comp.insert(el_no_comp.end(), nos_elem_vec[el_ref.nos[0]].begin(), nos_elem_vec[el_ref.nos[0]].end());
		el_no_comp.insert(el_no_comp.end(), nos_elem_vec[el_ref.nos[1]].begin(), nos_elem_vec[el_ref.nos[1]].end());
		el_no_comp.insert(el_no_comp.end(), nos_elem_vec[el_ref.nos[2]].begin(), nos_elem_vec[el_ref.nos[2]].end());

		sort(el_no_comp.begin(), el_no_comp.end());
		el_no_comp.erase(unique(el_no_comp.begin(), el_no_comp.end()), el_no_comp.end());
		
		bool no0 = false;
		bool no1 = false;
		bool no2 = false;		
		vector<Elem_Oposto> eo_vec;
		
		for (int i = 0; i < el_no_comp.size(); i++)
		{
			if (el_no_comp[i] == j) continue;
			
			Elemento el = malha.r_elem(el_no_comp[i]);

			if (checar_oposicao(el_ref, el, 0, 1) == true)
			{
				Elem_Oposto eo;
				eo.no = el_ref.nos[2];
				eo.elem = el_no_comp[i];
				eo_vec.push_back(eo);
				no2 = true;
			}
			if (checar_oposicao(el_ref, el, 0, 2) == true)
			{
				Elem_Oposto eo;
				eo.no = el_ref.nos[1];
				eo.elem = el_no_comp[i];
				eo_vec.push_back(eo);
				no1 = true;
			}
			if (checar_oposicao(el_ref, el, 1, 2) == true)
			{
				Elem_Oposto eo;
				eo.no = el_ref.nos[0];
				eo.elem = el_no_comp[i];
				eo_vec.push_back(eo);
				no0 = true;
			}
		}

		if (eo_vec.size() < 3)
		{
			if (no0 == false)
			{
				Elem_Oposto eo;
				eo.no = el_ref.nos[0];
				eo.elem = -1;
				eo_vec.push_back(eo);
			}
			if (no1 == false)
			{
				Elem_Oposto eo;		
				eo.no = el_ref.nos[1];
				eo.elem = -1;
				eo_vec.push_back(eo);
			}
			if (no2 == false)
			{
				Elem_Oposto eo;
				eo.no = el_ref.nos[2];
				eo.elem = -1;
				eo_vec.push_back(eo);
			}
		}
		elem_opostos.push_back(eo_vec);
	}
}

bool Semi_Lagrangiano::checar_contorno(Malha& malha, Elemento elem, int a, int b)
{
	vector<unsigned long> no_contorno_vec = malha.r_no_contorno();

	bool no1 = false;
	bool no2 = false;
	for (unsigned int i = 0; i < no_contorno_vec.size(); i++)
	{
		if (elem.nos[a] == no_contorno_vec[i]) no1 = true;
		if (elem.nos[b] == no_contorno_vec[i]) no2 = true;
	}

	if (no1 == true && no2 == true) return true;
	else return false;
}

void Semi_Lagrangiano::mostrar_Auxiliar()
{
	for (unsigned int i = 0; i < nos_elem_vec.size(); i++)
	{
		cout << "No " << i << " faz parte dos seguintes elementos: ";
		for (unsigned int j = 0; j < nos_elem_vec[i].size(); j++)
		{
			cout << nos_elem_vec[i][j] << " ";
		}
		cout << endl;
	}

	for (unsigned int i = 0; i < elem_opostos.size(); i++)
	{
		cout << "O elemento " << i << " tem os seguintes opostos:" << endl;
		for (unsigned int j = 0; j < elem_opostos[i].size(); j++)
		{
			cout << "elem " << elem_opostos[i][j].elem << " oposto ao no " << elem_opostos[i][j].no << endl;
		}
		cout << endl;
	}
}

bool Semi_Lagrangiano::checar_oposicao(Elemento elem_ref, Elemento elem, int a, int b)
{
	bool no1 = false;
	bool no2 = false;
	for (int i = 0; i < 3; i++)
	{
		if (elem.nos[i] == elem_ref.nos[a]) no1 = true;
		if (elem.nos[i] == elem_ref.nos[b]) no2 = true;
	}

	if (no1 == true && no2 == true) return true;
	else return false;
}

void Semi_Lagrangiano::semi_lagrangiano(Malha& malha, VectorXd& vx, VectorXd& vy, 
	double delta_t, VectorXd& prop, bool apenas_vert)
{
	system_clock::time_point total1 = system_clock::now();
	contagem_busca_linear = 0;
	duracao_busca_linear = 0;

	unsigned long gl = malha.r_num_nos();
	unsigned long vert = malha.r_num_vertices();
	unsigned long lim = gl;
	if (apenas_vert == true)  lim = vert;
	VectorXd prop_novo(gl);

	for (unsigned int i = 0; i < gl; i++)
	{		
		No no = malha.r_no(i);
		double x = no.x - vx[i] * delta_t;
		double y = no.y - vy[i] * delta_t;

		//cout << "--- NO " << i + 1 << " ---" << endl;
		//cout << "A posicao inical e: " << no.x << " " << no.y << endl;
		//cout << "A posicao e: " << x << " " << y << endl;
		//cout << "O vetor velocidade e: " << vx[i] << " " << vy[i] << endl;
		
		prop_novo[i] = achar_novo_valor(i, malha, x, y, prop);

		//cout << i << " " << prop_novo[i] << endl;
		//cout << endl;
		//system("pause");


		//long elem_index = busca_linear(malha, x, y);
		//double valor = 0;
		//if (elem_index == -1)
		//{
		//	unsigned long no_prox = buscar_no_prox2(malha, x, y);
		//	valor = prop[no_prox];
		//}
		//else valor = interpolar_Tri_Linear(malha, x, y, elem_index, prop);
		//prop_novo[i] = valor;

	}

	if (apenas_vert == true)
	{
		for (unsigned long i = 0; i < malha.r_num_elem(); i++)
		{
			Elemento el = malha.r_elem(i);
			unsigned long no3_index = el.nos[3]; //0 a 3, 4 nós
			No no3 = malha.r_no(no3_index);
			prop_novo[no3_index] = interpolar_Tri_Linear(malha, no3.x, no3.y, i, prop_novo);
		}
	}
	
	prop = prop_novo;

	cout << "O tempo das buscas lineares foi de " << duracao_busca_linear << endl;
	system_clock::time_point total2 = system_clock::now();
	cout << "O tempo total do SL foi de " << duration_cast<milliseconds>(total2 - total1).count() << " milisegundos" << endl;

	//cout << "A taxa percentual de tempo de busca linear por tempo total foi: " 
	//	<< (duracao_busca_linear / duration_cast<milliseconds>(total2 - total1).count()) * 100 << "%" << endl;
}

double Semi_Lagrangiano::achar_novo_valor(unsigned long index, Malha& malha, double x, double y, const VectorXd& prop)
{
	bool debug = false;
	long offset = 31;

	if (debug == true) cout << "Ponto: " << index+1 << endl;
	//if (index == 13) system("pause");

	vector<unsigned long> elems_do_no = nos_elem_vec[index];
	if (elems_do_no.empty() == true) cout << "Erro na funcao buscar elemento, vector vazio" << endl;

	//Escolhemos um dos elementos ao qual o nó pertence (qualquer um serve)
	long elem_index = elems_do_no[0];
	long elem_index_anterior = -1;

	if (debug == true) cout << "O primeiro elemento escolhido foi " << elem_index+offset << endl;

	//Checamos se o nó caiu dentro do primeiro elemento que escolhemos
	bool dentro = checar_ponto_dentro_elem(elem_index, x, y, malha);

	if (dentro == true)
	{
		if (debug == true) cout << "O ponto estava dentro do primeiro elemento testado" << endl;
		return interpolar_Tri_Linear(malha, x, y, elem_index, prop);
	}
	
	if (debug == true) cout << "O ponto --- NAO --- estava dentro do primeiro elemento testado, iniciando busca" << endl;

	//Se não, procuramos o próximo elemento
	while (true)
	{
		//Pegamos a lista de elementos opostos ao elemento escolhido
		vector<Elem_Oposto> elem_op = elem_opostos[elem_index];

		//Aqui checamos qual elemento da lista de elementos
		//opostos é o mais distante do nó de partida
		No no = malha.r_no(elem_op[0].no);
		long elem_oposto = elem_op[0].elem;
		double dist_ref = pow(x - no.x, 2) + pow(y - no.y, 2);

		for (unsigned int i = 1; i < elem_op.size(); i++)
		{
			No no = malha.r_no(elem_op[i].no);
			double dist = pow(x - no.x, 2) + pow(y - no.y, 2);
			if (dist > dist_ref)
			{
				dist_ref = dist;
				elem_oposto = elem_op[i].elem;
			}
		}

		if (elem_oposto == -1)
		{			
			No no_a = malha.r_no(index);
			double valor = interpolar_contorno(malha, elem_index, no_a.x, no_a.y, x, y, prop);
			if (debug == true) cout << "A interpolacao no contorno retornou " << valor << endl;
			return valor;
		}
		if (checar_ponto_dentro_elem(elem_oposto, x, y, malha) == true)
		{
			if (debug == true) cout << "Retornando elemento de index " << elem_oposto+offset << endl;
			return interpolar_Tri_Linear(malha, x, y, elem_oposto, prop);
		}
		if (elem_oposto == elem_index_anterior)
		{
			double novo_index = busca_linear(malha, x, y);
			if (debug == true) cout << "Busca linear retornou elemento " << elem_index+offset << endl;
			if (novo_index != -1) return interpolar_Tri_Linear(malha, x, y, novo_index, prop);
			else return interpolar_Tri_Linear(malha, x, y, elem_index, prop);
		}
		else
		{
			elem_index_anterior = elem_index;
			elem_index = elem_oposto;
			if (debug == true) cout << "novo index " << elem_index+ offset << endl;
		}
	}
}

double Semi_Lagrangiano::interpolar_contorno(Malha& malha, unsigned long elem_index, double xo, double yo, double x, double y, const VectorXd& prop)
{	
	vector<Elem_Oposto> elem_op = elem_opostos[elem_index];
	if (elem_op.size() != 3) cout << "PROBLEMA INTERPOLAR CONTORNO, ELEM_OP.SIZE != 3" << endl;

	vector<long> no_index_vec;
	for (int i = 0; i < elem_op.size(); i++) if (elem_op[i].elem != -1) no_index_vec.push_back(elem_op[i].no);
	if (no_index_vec.size() != 2) cout << "PROBLEMA INTERPOLAR CONTORNO, no_index_vec.size() != 2" << endl;

	No a = malha.r_no(no_index_vec[0]);
	No b = malha.r_no(no_index_vec[1]);
	
	double d = (xo - x) * (a.y - b.y) - (yo - y) * (a.x - b.x);

	//if (d == 0)	cout << "PROBLEMA INTERPOLAR, D = 0" << endl;
	double px = ((xo * y - yo * x) * (a.x - b.x) - (xo - x) * (a.x * b.y - a.y * b.x)) / d;
	double py = ((xo * y - yo * x) * (a.y - b.y) - (yo - y) * (a.x * b.y - a.y * b.x)) / d;

	double Lpb = sqrt(pow(px - b.x, 2) + pow(py - b.y, 2));
	double Lab = sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));

	double peso_a = Lpb / Lab;
	double peso_b = 1 - peso_a;

	if (peso_a < 0) peso_a = 0;
	if (peso_b < 0) peso_b = 0;
	if (peso_a > 1) peso_a = 1;
	if (peso_b > 1) peso_b = 1;

	if (peso_a > 1 || peso_a < 0 || peso_b > 1 || peso_b < 0)
	{
		cout << "PROBLEMA INTERPOLAR CONTORNO, pesos incoerentes" << endl;
	}

	double valor = prop[a.id] * peso_a + prop[b.id] * peso_b;
	if (isnan(valor) || d == 0) valor = prop[buscar_no_prox2(malha, x, y)];
	return valor;
	
	
	////Esta função pode ser melhorada depois. Por hora, vou só buscar o ponto mais
	////próximo, mas tem como achar o elemento, achar a intersecção e fazer a inter-
	////polação apropriada

	//Elemento el = malha.r_elem(elem_index);	
	//No no = malha.r_no(el.nos[0]);
	//unsigned long no_mais_prox = no.id;

	//double dist, dist_ref;
	//dist_ref = pow(x - no.x, 2) + pow(y - no.y, 2);

	//for (unsigned int i = 1; i < el.nos.size(); i++)
	//{
	//	No no = malha.r_no(el.nos[i]);
	//	//Como é apenas pra comparação, e não há interesse em achar o valor
	//	//real da distância, pulei o cálculo da raiz, que tem custo computacional
	//	dist = pow(x - no.x, 2) + pow(y - no.y, 2);

	//	if (dist < dist_ref)
	//	{
	//		dist_ref = dist;
	//		no_mais_prox = no.id;
	//	}
	//}
	//return no_mais_prox;
}

long Semi_Lagrangiano::busca_linear(Malha& malha, double x, double y)
{
	//Busca linear, percorrendo todos elementos
	for (unsigned int i = 0; i < malha.r_num_elem(); i++)
	{
		if (checar_ponto_dentro_elem(i, x, y, malha) == true) return i;
	}
	return -1;
}

double Semi_Lagrangiano::area_triangulo(double x1, double y1, double x2, double y2, double x3, double y3)
{
	Matrix3d mat_A;
	mat_A << 1, x1, y1, 1, x2, y2, 1, x3, y3;

	double A = mat_A.determinant() * 0.5;
	return A;
}

double Semi_Lagrangiano::interpolar_Tri_Linear(Malha& malha, double x, double y,
	long elem_index, const VectorXd& prop)
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

bool Semi_Lagrangiano::checar_ponto_dentro_elem(unsigned long elem_index, double x, double y, Malha& malha)
{
	Elemento elem = malha.r_elem(elem_index);

	No no1 = malha.r_no(elem.nos[0]);
	No no2 = malha.r_no(elem.nos[1]);
	No no3 = malha.r_no(elem.nos[2]);

	double A = area_triangulo(no1.x, no1.y, no2.x, no2.y, no3.x, no3.y);
	double A1 = area_triangulo(x, y, no2.x, no2.y, no3.x, no3.y);
	double A2 = area_triangulo(no1.x, no1.y, x, y, no3.x, no3.y);
	double A3 = area_triangulo(no1.x, no1.y, no2.x, no2.y, x, y);

	if (A <= 0 && A1 <= 0 && A2 <= 0 && A3 <= 0) return true;
	else if ((A >= 0 && A1 >= 0 && A2 >= 0 && A3 >= 0)) return true;
	else return false;
	
	//if (A1 > A) return false;
	//if (A2 > A) return false;
	//if (A3 > A) return false;
	//return true;

	//double soma = abs(A1) + abs(A2) + abs(A3);
	//double percentual = soma / A;

	//if (percentual >= 0.999999 && percentual <= 1.00001) return true;
	//else return false;
}

//TEMPORARIO para testes
unsigned long Semi_Lagrangiano::buscar_no_prox2(Malha& malha, double x, double y)
{
	//Esta função pode ser melhorada depois. Por hora, vou só buscar o ponto mais
	//próximo, mas tem como achar o elemento, achar a intersecção e fazer a inter-
	//polação apropriada

	vector<unsigned long> nos_contorno = malha.r_no_contorno();

	double dist, dist_ref;
	unsigned long no_prox = nos_contorno[0];

	No no = malha.r_no(nos_contorno[0]);
	dist_ref = pow(x - no.x, 2) + pow(y - no.y, 2);

	for (unsigned int i = 1; i < nos_contorno.size(); i++)
	{
		No no = malha.r_no(nos_contorno[i]);

		//Como é apenas pra comparação, e não há interesse em achar o valor
		//real da distância, pulei o cálculo da raiz, que tem custo computacional
		dist = pow(x - no.x, 2) + pow(y - no.y, 2);

		if (dist < dist_ref)
		{
			dist_ref = dist;
			no_prox = nos_contorno[i];
		}
	}
	return no_prox;
}