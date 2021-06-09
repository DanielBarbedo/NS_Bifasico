#include "NavierStokes2D.h"

NavierStokes2D::NavierStokes2D(string nome)
{
	malha.ler_malha(nome);	
	malha.gerar_elemento_mini();
	//malha.mostrar_malha();

	gl = malha.r_num_nos();
	vert = malha.r_num_vertices();

	cout << "--- Malha lida, com " << gl << " nos e " << malha.r_num_elem() << " elementos. ---" << endl << endl;
	cout << "Alocando espaco para matrizes" << endl << endl;
	
	K = SparseMatrix<double>(gl, gl);
	M = SparseMatrix<double>(gl, gl);
	Gx = SparseMatrix<double>(gl, vert);
	Gy = SparseMatrix<double>(gl, vert);

	//---Avaliacao de propriedade---
	gerar_linha_prop(0, 0, 0, 0, 0);
	//---Geração de arquivos de saída
	max_saida = 400;
	iteracoes_por_saida = 1;
}

void NavierStokes2D::resolver_stokes_permanente()
{	
	cout << "Montando matrizes fixas" << endl;	
	system_clock::time_point t1 = system_clock::now();	
	montar_matrizes_fixas(K, M, Gx, Gy);
	system_clock::time_point t2 = system_clock::now();
	cout << "A montagem das matrizes fixas levou " << duration_cast<seconds>(t2 - t1).count() << " segundos" << endl << endl;

	cout << "Montando matriz A a partir das matrizes fixas" << endl;
	t1 = system_clock::now();
	VectorXd B(gl * 2 + vert);
	VectorXd sol(gl * 2 + vert);
	B.fill(0);
	SparseMatrix<double, Eigen::RowMajor> A(2 * gl + vert, 2 * gl + vert);	
	vector<Triplet<double>> A_triplets = montar_matriz_A(false, 1, 1);
	A.setFromTriplets(A_triplets.begin(), A_triplets.end());
	t2 = system_clock::now();
	cout << "O assembly de A levou " << duration_cast<seconds>(t2 - t1).count() << " segundos" << endl << endl;

	cout << "Aplicando condicoes de contorno a matriz A e ao vetor B" << endl;
	t1 = system_clock::now();
	aplicar_cc_A(A);
	aplicar_cc_B(B);
	t2 = system_clock::now();
	cout << "A aplicacao das condicoes de contorno a matriz A e ao vetor B levou " << duration_cast<seconds>(t2 - t1).count() << " segundos" << endl << endl;

	cout << "Resolvendo sistema linear" << endl;
	t1 = system_clock::now();
	//Eigen::BiCGSTAB<SparseMatrix<double>> solver;
	Eigen::SparseLU<SparseMatrix<double>> solver;
	solver.compute(A);
	if (solver.info() != EXIT_SUCCESS) cout << "Decomposicao falhou" << endl;
	sol = solver.solve(B);
	t2 = system_clock::now();
	cout << "A resolucao do sistema linear levou " << duration_cast<seconds>(t2 - t1).count() << " segundos" << endl << endl;

	cout << "Gerando arquivo de saída" << endl;
	t1 = system_clock::now();
	gerar_arquivo_saida(sol, 0, 1, 1, "Perm");
	t2 = system_clock::now();
	cout << "A geracao do arquivo de saida levou " << duration_cast<seconds>(t2 - t1).count() << " segundos" << endl << endl;
}

void NavierStokes2D::resolver_stokes_transiente(double delta_t, unsigned long max_iter)
{
	iteracoes_por_saida = max_iter / max_saida;
	if (iteracoes_por_saida < 1) iteracoes_por_saida = 1;

	cout << "Montando matrizes fixas" << endl;	
	montar_matrizes_fixas(K, M, Gx, Gy);

	cout << "Montando matriz A a partir das matrizes fixas" << endl;
	SparseMatrix<double, Eigen::RowMajor> A(2 * gl + vert, 2 * gl + vert);
	vector<Triplet<double>> A_triplets = montar_matriz_A(true, delta_t, 1);
	A.setFromTriplets(A_triplets.begin(), A_triplets.end());
	
	cout << "Aplicando condicoes de contorno a matriz A" << endl;
	aplicar_cc_A(A);

	cout << "Iniciando processo iterativo" << endl;
	VectorXd B(gl * 2 + vert);
	VectorXd sol(gl * 2 + vert);
	VectorXd M_vx(gl);
	VectorXd M_vy(gl);
	M_vx.fill(0);
	M_vy.fill(0);
	for (unsigned int i = 0; i < max_iter; i++)
	{
		cout << "Iteracao: " << i << endl;
		B.fill(0);
		M_vx = M * ((double)1 / delta_t) * sol.head(gl);
		M_vy = M * ((double)1 / delta_t) * sol.segment(gl, gl);
		B.head(gl) = M_vx;
		B.segment(gl, gl) = M_vy;

		cout << "Aplicando condicoes de contorno ao vetor B" << endl;
		aplicar_cc_B(B);

		cout << "Resolvendo sistema linear" << endl;
		//Eigen::BiCGSTAB<SparseMatrix<double>> solver;
		Eigen::SparseLU<SparseMatrix<double>> solver;
		solver.compute(A);
		if (solver.info() != EXIT_SUCCESS) cout << "Decomposicao falhou" << endl;
		sol = solver.solve(B);

		gerar_arquivo_saida(sol, i, delta_t, 1, "Transiente");		
	}	
}

void NavierStokes2D::resolver_navier_stokes(double delta_t, unsigned long max_iter, double Re, string tag)
{
	iteracoes_por_saida = max_iter / max_saida;
	if (iteracoes_por_saida < 1) iteracoes_por_saida = 1;
	
	system_clock::time_point total1 = system_clock::now();
	
	cout << "Montando matrizes fixas" << endl;
	system_clock::time_point t1 = system_clock::now();
	montar_matrizes_fixas(K, M, Gx, Gy);	
	system_clock::time_point t2 = system_clock::now();
	cout << "A montagem das matrizes fixas levou " << duration_cast<seconds>(t2 - t1).count() << " segundos" << endl << endl;

	cout << "Montando A a partir das matrizes fixas" << endl;
	t1 = system_clock::now();
	SparseMatrix<double, Eigen::RowMajor> A(2 * gl + vert, 2 * gl + vert);
	vector<Triplet<double>> A_triplets = montar_matriz_A(true, delta_t, Re);
	A.setFromTriplets(A_triplets.begin(), A_triplets.end());
	t2 = system_clock::now();
	cout << "A montagem de A levou " << duration_cast<seconds>(t2 - t1).count() << " segundos" << endl << endl;

	cout << "Aplicando ccs na matriz A" << endl;
	t1 = system_clock::now();
	aplicar_cc_A(A);
	t2 = system_clock::now();
	cout << "A aplicação das ccs em A levou " << duration_cast<seconds>(t2 - t1).count() << " segundos" << endl << endl;

	cout << "Iniciando o semi Lagrangiano" << endl;
	t1 = system_clock::now();
	Semi_Lagrangiano SL(malha);
	t2 = system_clock::now();
	cout << "A inicializacao do semi Lagrangiano levou " << duration_cast<seconds>(t2 - t1).count() << " segundos" << endl << endl;

	cout << "Iniciando processo iterativo" << endl << endl;
	VectorXd B(gl * 2 + vert);
	VectorXd sol(gl * 2 + vert);	
	VectorXd vx(gl);
	VectorXd vy(gl);
	//VectorXd vx_temp(gl);
	//VectorXd vy_temp(gl);
	//vx.fill(0);
	//vy.fill(0);
	sol.fill(0);

	for (unsigned int i = 0; i < max_iter; i++)
	{
		cout << "Iteracao: " << i << endl;
		t1 = system_clock::now();	

		vx = sol.head(gl);
		vy = sol.segment(gl, gl);
		cout << "Utilizando semi lagrangiano" << endl;		
		SL.semi_lagrangiano(malha, vx, vy, delta_t);
		SL.multiplicar_SL(vx); //Plotar as velocidades convectadas
		SL.multiplicar_SL(vy);

		B.fill(0);
		B.head(gl) = M * ((double)1 / delta_t) * vx;
		B.segment(gl, gl) = M * ((double)1 / delta_t) * vy;

		cout << "Aplicando condicoes de contorno ao vetor B" << endl;
		aplicar_cc_B(B);

		cout << "Resolvendo sistema linear" << endl << endl;
		//Eigen::BiCGSTAB<SparseMatrix<double>> solver;
		Eigen::SparseLU<SparseMatrix<double>> solver;
		solver.compute(A);
		if (solver.info() != EXIT_SUCCESS) cout << "Decomposicao falhou" << endl;
		sol = solver.solve(B);

		gerar_arquivo_saida(sol, i, delta_t, Re, tag);
		t2 = system_clock::now();
		cout << "A iteracao " << i << " levou " << duration_cast<seconds>(t2 - t1).count() << " segundos" << endl << endl;

		//system("pause");
	}

	avaliar_prop(sol, SL, tag, linha_prop_num_pontos, linha_prop_xo, linha_prop_yo, linha_prop_x, linha_prop_y);

	system_clock::time_point total2 = system_clock::now();
	cout << "A malha possui " << gl << " nos, e " << malha.r_num_elem() << " elementos." << endl;
	cout << "Todo procedimento (excluindo a leitura da malha) levou " << duration_cast<seconds>(total2 - total1).count() << " segundos" << endl;
}

//----------

void NavierStokes2D::montar_matrizes_fixas(SparseMatrix<double>& K, SparseMatrix<double>& M, SparseMatrix<double>& Gx, SparseMatrix<double>& Gy)
{
	vector<Triplet<double>> K_triplets;
	vector<Triplet<double>> M_triplets;
	vector<Triplet<double>> Gx_triplets;
	vector<Triplet<double>> Gy_triplets;

	for (unsigned int i = 0; i < malha.r_num_elem(); i++)
	{
		Elemento elem = malha.r_elem(i);

		No no1 = malha.r_no(elem.nos[0]);
		No no2 = malha.r_no(elem.nos[1]);
		No no3 = malha.r_no(elem.nos[2]);
		No no4 = malha.r_no(elem.nos[3]);

		//Contas Auxiliares
		double ai = no2.x * no3.y - no3.x * no2.y;
		double aj = no3.x * no1.y - no1.x * no3.y;
		double ak = no1.x * no2.y - no2.x * no1.y;

		double bi = no2.y - no3.y;
		double bj = no3.y - no1.y;
		double bk = no1.y - no2.y;

		double ci = no3.x - no2.x;
		double cj = no1.x - no3.x;
		double ck = no2.x - no1.x;

		Matrix3d mat_A;
		mat_A << 1, no1.x, no1.y, 1, no2.x, no2.y, 1, no3.x, no3.y;
		double A = mat_A.determinant() * 0.5;

		double z = (bj * bj + bj * bk + bk * bk + cj * cj + cj * ck + ck * ck) / (4 * A);

		Matrix3d mat_kxe_lin;
		mat_kxe_lin(0, 0) = bi * bi;
		mat_kxe_lin(0, 1) = bi * bj;
		mat_kxe_lin(0, 2) = bi * bk;
		mat_kxe_lin(1, 0) = bj * bi;
		mat_kxe_lin(1, 1) = bj * bj;
		mat_kxe_lin(1, 2) = bj * bk;
		mat_kxe_lin(2, 0) = bk * bi;
		mat_kxe_lin(2, 1) = bk * bj;
		mat_kxe_lin(2, 2) = bk * bk;
		mat_kxe_lin = mat_kxe_lin / (4 * A);

		Matrix3d mat_kye_lin;
		mat_kye_lin(0, 0) = ci * ci;
		mat_kye_lin(0, 1) = ci * cj;
		mat_kye_lin(0, 2) = ci * ck;
		mat_kye_lin(1, 0) = cj * ci;
		mat_kye_lin(1, 1) = cj * cj;
		mat_kye_lin(1, 2) = cj * ck;
		mat_kye_lin(2, 0) = ck * ci;
		mat_kye_lin(2, 1) = ck * cj;
		mat_kye_lin(2, 2) = ck * ck;
		mat_kye_lin = mat_kye_lin / (4 * A);

		Matrix3d mat_ke_lin;
		mat_ke_lin = mat_kxe_lin + mat_kye_lin;

		MatrixXd mat_ke_mini(4, 4);
		mat_ke_mini(0, 0) = mat_ke_lin(0, 0) + z * 9 / 10;
		mat_ke_mini(0, 1) = mat_ke_lin(0, 1) + z * 9 / 10;
		mat_ke_mini(0, 2) = mat_ke_lin(0, 2) + z * 9 / 10;
		mat_ke_mini(0, 3) = z * -27 / 10;
		mat_ke_mini(1, 0) = mat_ke_lin(1, 0) + z * 9 / 10;
		mat_ke_mini(1, 1) = mat_ke_lin(1, 1) + z * 9 / 10;
		mat_ke_mini(1, 2) = mat_ke_lin(1, 2) + z * 9 / 10;
		mat_ke_mini(1, 3) = z * -27 / 10;
		mat_ke_mini(2, 0) = mat_ke_lin(2, 0) + z * 9 / 10;
		mat_ke_mini(2, 1) = mat_ke_lin(2, 1) + z * 9 / 10;
		mat_ke_mini(2, 2) = mat_ke_lin(2, 2) + z * 9 / 10;
		mat_ke_mini(2, 3) = z * -27 / 10;
		mat_ke_mini(3, 0) = z * -27 / 10;
		mat_ke_mini(3, 1) = z * -27 / 10;
		mat_ke_mini(3, 2) = z * -27 / 10;
		mat_ke_mini(3, 3) = z * 81 / 10;

		//cout << "Mat_ke_lin:" << endl << mat_ke_lin << endl << endl;
		//cout << "Mat_ke_mini:" << endl << mat_ke_mini << endl << endl;

		MatrixXd mat_me_mini(4, 4);
		mat_me_mini(0, 0) = 83;
		mat_me_mini(0, 1) = 13;
		mat_me_mini(0, 2) = 13;
		mat_me_mini(0, 3) = 45;
		mat_me_mini(1, 0) = 13;
		mat_me_mini(1, 1) = 83;
		mat_me_mini(1, 2) = 13;
		mat_me_mini(1, 3) = 45;
		mat_me_mini(2, 0) = 13;
		mat_me_mini(2, 1) = 13;
		mat_me_mini(2, 2) = 83;
		mat_me_mini(2, 3) = 45;
		mat_me_mini(3, 0) = 45;
		mat_me_mini(3, 1) = 45;
		mat_me_mini(3, 2) = 45;
		mat_me_mini(3, 3) = 243;
		mat_me_mini = mat_me_mini * A / 840;

		Matrix3d mat_gxe_lin;
		mat_gxe_lin << bi, bj, bk, bi, bj, bk, bi, bj, bk;
		mat_gxe_lin = mat_gxe_lin / 6;
		Matrix3d mat_gxe_lin_trans;
		mat_gxe_lin_trans = mat_gxe_lin.transpose();

		Matrix3d mat_gye_lin;
		mat_gye_lin << ci, cj, ck, ci, cj, ck, ci, cj, ck;
		mat_gye_lin = mat_gye_lin / 6;
		Matrix3d mat_gye_lin_trans;
		mat_gye_lin_trans = mat_gye_lin.transpose();

		MatrixXd mat_gxe_mini(4, 3);
		mat_gxe_mini(0, 0) = mat_gxe_lin(0, 0) * ((double)-9 / 20) - mat_gxe_lin_trans(0, 0);
		mat_gxe_mini(0, 1) = mat_gxe_lin(0, 1) * ((double)-9 / 20) - mat_gxe_lin_trans(0, 1);
		mat_gxe_mini(0, 2) = mat_gxe_lin(0, 2) * ((double)-9 / 20) - mat_gxe_lin_trans(0, 2);
		mat_gxe_mini(1, 0) = mat_gxe_lin(1, 0) * ((double)-9 / 20) - mat_gxe_lin_trans(1, 0);
		mat_gxe_mini(1, 1) = mat_gxe_lin(1, 1) * ((double)-9 / 20) - mat_gxe_lin_trans(1, 1);
		mat_gxe_mini(1, 2) = mat_gxe_lin(1, 2) * ((double)-9 / 20) - mat_gxe_lin_trans(1, 2);
		mat_gxe_mini(2, 0) = mat_gxe_lin(2, 0) * ((double)-9 / 20) - mat_gxe_lin_trans(2, 0);
		mat_gxe_mini(2, 1) = mat_gxe_lin(2, 1) * ((double)-9 / 20) - mat_gxe_lin_trans(2, 1);
		mat_gxe_mini(2, 2) = mat_gxe_lin(2, 2) * ((double)-9 / 20) - mat_gxe_lin_trans(2, 2);
		mat_gxe_mini(3, 0) = ((double)9 / 40) * bi;
		mat_gxe_mini(3, 1) = ((double)9 / 40) * bj;
		mat_gxe_mini(3, 2) = ((double)9 / 40) * bk;

		MatrixXd mat_gye_mini(4, 3);
		mat_gye_mini(0, 0) = mat_gye_lin(0, 0) * ((double)-9 / 20) - mat_gye_lin_trans(0, 0);
		mat_gye_mini(0, 1) = mat_gye_lin(0, 1) * ((double)-9 / 20) - mat_gye_lin_trans(0, 1);
		mat_gye_mini(0, 2) = mat_gye_lin(0, 2) * ((double)-9 / 20) - mat_gye_lin_trans(0, 2);
		mat_gye_mini(1, 0) = mat_gye_lin(1, 0) * ((double)-9 / 20) - mat_gye_lin_trans(1, 0);
		mat_gye_mini(1, 1) = mat_gye_lin(1, 1) * ((double)-9 / 20) - mat_gye_lin_trans(1, 1);
		mat_gye_mini(1, 2) = mat_gye_lin(1, 2) * ((double)-9 / 20) - mat_gye_lin_trans(1, 2);
		mat_gye_mini(2, 0) = mat_gye_lin(2, 0) * ((double)-9 / 20) - mat_gye_lin_trans(2, 0);
		mat_gye_mini(2, 1) = mat_gye_lin(2, 1) * ((double)-9 / 20) - mat_gye_lin_trans(2, 1);
		mat_gye_mini(2, 2) = mat_gye_lin(2, 2) * ((double)-9 / 20) - mat_gye_lin_trans(2, 2);
		mat_gye_mini(3, 0) = ((double)9 / 40) * ci;
		mat_gye_mini(3, 1) = ((double)9 / 40) * cj;
		mat_gye_mini(3, 2) = ((double)9 / 40) * ck;

		//Incluindo os resultados locais nas matrizes globais
		K_triplets.push_back(Triplet<double>(no1.id, no1.id, mat_ke_mini(0, 0)));
		K_triplets.push_back(Triplet<double>(no1.id, no2.id, mat_ke_mini(0, 1)));
		K_triplets.push_back(Triplet<double>(no1.id, no3.id, mat_ke_mini(0, 2)));
		K_triplets.push_back(Triplet<double>(no1.id, no4.id, mat_ke_mini(0, 3)));
		K_triplets.push_back(Triplet<double>(no2.id, no1.id, mat_ke_mini(1, 0)));
		K_triplets.push_back(Triplet<double>(no2.id, no2.id, mat_ke_mini(1, 1)));
		K_triplets.push_back(Triplet<double>(no2.id, no3.id, mat_ke_mini(1, 2)));
		K_triplets.push_back(Triplet<double>(no2.id, no4.id, mat_ke_mini(1, 3)));
		K_triplets.push_back(Triplet<double>(no3.id, no1.id, mat_ke_mini(2, 0)));
		K_triplets.push_back(Triplet<double>(no3.id, no2.id, mat_ke_mini(2, 1)));
		K_triplets.push_back(Triplet<double>(no3.id, no3.id, mat_ke_mini(2, 2)));
		K_triplets.push_back(Triplet<double>(no3.id, no4.id, mat_ke_mini(2, 3)));
		K_triplets.push_back(Triplet<double>(no4.id, no1.id, mat_ke_mini(3, 0)));
		K_triplets.push_back(Triplet<double>(no4.id, no2.id, mat_ke_mini(3, 1)));
		K_triplets.push_back(Triplet<double>(no4.id, no3.id, mat_ke_mini(3, 2)));
		K_triplets.push_back(Triplet<double>(no4.id, no4.id, mat_ke_mini(3, 3)));

		M_triplets.push_back(Triplet<double>(no1.id, no1.id, mat_me_mini(0, 0)));
		M_triplets.push_back(Triplet<double>(no1.id, no2.id, mat_me_mini(0, 1)));
		M_triplets.push_back(Triplet<double>(no1.id, no3.id, mat_me_mini(0, 2)));
		M_triplets.push_back(Triplet<double>(no1.id, no4.id, mat_me_mini(0, 3)));
		M_triplets.push_back(Triplet<double>(no2.id, no1.id, mat_me_mini(1, 0)));
		M_triplets.push_back(Triplet<double>(no2.id, no2.id, mat_me_mini(1, 1)));
		M_triplets.push_back(Triplet<double>(no2.id, no3.id, mat_me_mini(1, 2)));
		M_triplets.push_back(Triplet<double>(no2.id, no4.id, mat_me_mini(1, 3)));
		M_triplets.push_back(Triplet<double>(no3.id, no1.id, mat_me_mini(2, 0)));
		M_triplets.push_back(Triplet<double>(no3.id, no2.id, mat_me_mini(2, 1)));
		M_triplets.push_back(Triplet<double>(no3.id, no3.id, mat_me_mini(2, 2)));
		M_triplets.push_back(Triplet<double>(no3.id, no4.id, mat_me_mini(2, 3)));
		M_triplets.push_back(Triplet<double>(no4.id, no1.id, mat_me_mini(3, 0)));
		M_triplets.push_back(Triplet<double>(no4.id, no2.id, mat_me_mini(3, 1)));
		M_triplets.push_back(Triplet<double>(no4.id, no3.id, mat_me_mini(3, 2)));
		M_triplets.push_back(Triplet<double>(no4.id, no4.id, mat_me_mini(3, 3)));

		Gx_triplets.push_back(Triplet<double>(no1.id, no1.id, mat_gxe_mini(0, 0)));
		Gx_triplets.push_back(Triplet<double>(no1.id, no2.id, mat_gxe_mini(0, 1)));
		Gx_triplets.push_back(Triplet<double>(no1.id, no3.id, mat_gxe_mini(0, 2)));
		Gx_triplets.push_back(Triplet<double>(no2.id, no1.id, mat_gxe_mini(1, 0)));
		Gx_triplets.push_back(Triplet<double>(no2.id, no2.id, mat_gxe_mini(1, 1)));
		Gx_triplets.push_back(Triplet<double>(no2.id, no3.id, mat_gxe_mini(1, 2)));
		Gx_triplets.push_back(Triplet<double>(no3.id, no1.id, mat_gxe_mini(2, 0)));
		Gx_triplets.push_back(Triplet<double>(no3.id, no2.id, mat_gxe_mini(2, 1)));
		Gx_triplets.push_back(Triplet<double>(no3.id, no3.id, mat_gxe_mini(2, 2)));
		Gx_triplets.push_back(Triplet<double>(no4.id, no1.id, mat_gxe_mini(3, 0)));
		Gx_triplets.push_back(Triplet<double>(no4.id, no2.id, mat_gxe_mini(3, 1)));
		Gx_triplets.push_back(Triplet<double>(no4.id, no3.id, mat_gxe_mini(3, 2)));

		Gy_triplets.push_back(Triplet<double>(no1.id, no1.id, mat_gye_mini(0, 0)));
		Gy_triplets.push_back(Triplet<double>(no1.id, no2.id, mat_gye_mini(0, 1)));
		Gy_triplets.push_back(Triplet<double>(no1.id, no3.id, mat_gye_mini(0, 2)));
		Gy_triplets.push_back(Triplet<double>(no2.id, no1.id, mat_gye_mini(1, 0)));
		Gy_triplets.push_back(Triplet<double>(no2.id, no2.id, mat_gye_mini(1, 1)));
		Gy_triplets.push_back(Triplet<double>(no2.id, no3.id, mat_gye_mini(1, 2)));
		Gy_triplets.push_back(Triplet<double>(no3.id, no1.id, mat_gye_mini(2, 0)));
		Gy_triplets.push_back(Triplet<double>(no3.id, no2.id, mat_gye_mini(2, 1)));
		Gy_triplets.push_back(Triplet<double>(no3.id, no3.id, mat_gye_mini(2, 2)));
		Gy_triplets.push_back(Triplet<double>(no4.id, no1.id, mat_gye_mini(3, 0)));
		Gy_triplets.push_back(Triplet<double>(no4.id, no2.id, mat_gye_mini(3, 1)));
		Gy_triplets.push_back(Triplet<double>(no4.id, no3.id, mat_gye_mini(3, 2)));
	}

	K.setFromTriplets(K_triplets.begin(), K_triplets.end());
	M.setFromTriplets(M_triplets.begin(), M_triplets.end());
	Gx.setFromTriplets(Gx_triplets.begin(), Gx_triplets.end());
	Gy.setFromTriplets(Gy_triplets.begin(), Gy_triplets.end());
}

vector<Triplet<double>> NavierStokes2D::montar_matriz_A(bool transiente, double delta_t, double Re)
{
	vector<Triplet<double>> A_triplets;
	for (int i = 0; i < K.outerSize(); ++i)
	{
		for (SparseMatrix<double>::InnerIterator it(K, i); it; ++it)
		{
			//it.value();
			//it.row();   // row index
			//it.col();   // col index (here it is equal to k)
			//it.index(); // inner index, here it is equal to it.row()
			A_triplets.push_back(Triplet<double>((long)it.row(), (long)it.col(), it.value() / Re));
			A_triplets.push_back(Triplet<double>((long)it.row() + gl, (long)it.col() + gl, it.value() / Re));
		}
	}
	for (int i = 0; i < Gx.outerSize(); ++i)
	{
		for (SparseMatrix<double>::InnerIterator it(Gx, i); it; ++it)
		{
			A_triplets.push_back(Triplet<double>((long)it.row(), (long)it.col() + gl * 2, it.value()));
			A_triplets.push_back(Triplet<double>((long)it.col() + gl * 2, (long)it.row(), -it.value()));
		}
	}
	for (int i = 0; i < Gy.outerSize(); ++i)
	{
		for (SparseMatrix<double>::InnerIterator it(Gy, i); it; ++it)
		{
			A_triplets.push_back(Triplet<double>((long)it.row() + gl, (long)it.col() + gl * 2, it.value()));
			A_triplets.push_back(Triplet<double>((long)it.col() + gl * 2, (long)it.row() + gl, -it.value()));
		}
	}
	if (transiente == true)
	{
		for (int i = 0; i < M.outerSize(); ++i)
		{
			for (SparseMatrix<double>::InnerIterator it(M, i); it; ++it)
			{
				A_triplets.push_back(Triplet<double>((long)it.row(), (long)it.col(), it.value() / delta_t));
				A_triplets.push_back(Triplet<double>((long)it.row() + gl, (long)it.col() + gl, it.value() / delta_t));
			}
		}
	}
	return A_triplets;
}

void NavierStokes2D::aplicar_cc_A(SparseMatrix<double, Eigen::RowMajor>& A)
{
	unordered_map<unsigned long, double> cc_vx_map = malha.r_cc_vx();
	for (unordered_map<unsigned long, double>::iterator it_map = cc_vx_map.begin(); it_map != cc_vx_map.end(); it_map++)
	{
		for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(A, it_map->first); it; ++it)
		{
			//it.valueRef() = 0;
			A.coeffRef(it.row(), it.col()) = 0;
		}
		A.coeffRef(it_map->first, it_map->first) = 1;
	}	
	unordered_map<unsigned long, double> cc_vy_map = malha.r_cc_vy();
	for (unordered_map<unsigned long, double>::iterator it_map = cc_vy_map.begin(); it_map != cc_vy_map.end(); it_map++)
	{
		for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(A, it_map->first + gl); it; ++it)
		{
			//it.valueRef() = 0;
			A.coeffRef(it.row(), it.col()) = 0;
		}
		A.coeffRef(it_map->first + gl, it_map->first + gl) = 1;
	}
	unordered_map<unsigned long, double> cc_p_map = malha.r_cc_p();
	for (unordered_map<unsigned long, double>::iterator it_map = cc_p_map.begin(); it_map != cc_p_map.end(); it_map++)
	{
		for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(A, it_map->first + gl * 2); it; ++it)
		{
			//it.valueRef() = 0;
			A.coeffRef(it.row(), it.col()) = 0;
		}
		A.coeffRef(it_map->first + gl * 2, it_map->first + gl * 2) = 1;
	}
}

void NavierStokes2D::aplicar_cc_B(VectorXd& B)
{
	unordered_map<unsigned long, double> cc_vx_map = malha.r_cc_vx();
	for (unordered_map<unsigned long, double>::iterator it = cc_vx_map.begin(); it != cc_vx_map.end(); it++)
	{
		B(it->first) = it->second;
	}
	unordered_map<unsigned long, double> cc_vy_map = malha.r_cc_vy();
	for (unordered_map<unsigned long, double>::iterator it = cc_vy_map.begin(); it != cc_vy_map.end(); it++)
	{
		B(it->first + gl) = it->second;
	}
	unordered_map<unsigned long, double> cc_p_map = malha.r_cc_p();
	for (unordered_map<unsigned long, double>::iterator it = cc_p_map.begin(); it != cc_p_map.end(); it++)
	{
		B(it->first + gl * 2) = it->second;
	}
}

void NavierStokes2D::gerar_arquivo_saida(VectorXd& sol, unsigned long iter, double delta_t, double Re, string tag)
{
	if (iter % iteracoes_por_saida != 0) return;
	
	fstream arquivo_saida;

	string Re_string = to_string(Re);
	if(Re >= 1) Re_string = Re_string.substr(0, Re_string.size() - 7);
	else Re_string = Re_string.substr(0, Re_string.size() - 2);
	string dt_string = to_string(delta_t);
	dt_string = dt_string.substr(0, dt_string.size() - 2);	
	
	arquivo_saida.open(malha.r_nome() + "_dt_" + dt_string +"_Re_" + Re_string
		+ "_" + tag + "_" + to_string(iter) + ".vtk", ios::out);
	unsigned long gl = malha.r_num_nos();
	unsigned long vert = malha.r_num_vertices();

	arquivo_saida << "# vtk DataFile Version 1.0" << endl;
	arquivo_saida << "Funcao Navier-Stokes 2D" << endl;
	arquivo_saida << "ASCII" << endl;
	arquivo_saida << "DATASET UNSTRUCTURED_GRID" << endl;
	arquivo_saida << endl;

	arquivo_saida << "POINTS " << vert << " double" << endl;
	for (unsigned int i = 0; i < vert; i++)
	{
		No no = malha.r_no(i);
		arquivo_saida << no.x << " " << no.y << " " << no.z << endl;
	}
	arquivo_saida << endl;

	unsigned long tamanho_lista = 0;
	for (unsigned long i = 0; i < malha.r_num_elem(); i++)
	{
		int num_nos = 0;
		Elemento el = malha.r_elem(i);
		if (el.tipo == Tipo_Elem::Tri_Linear) num_nos = 3;
		if (el.tipo == Tipo_Elem::Mini) num_nos = 3;
		if (el.tipo == Tipo_Elem::Tri_Quadr) num_nos = 6;
		tamanho_lista += 1 + num_nos;
	}
	arquivo_saida << "CELLS " << malha.r_num_elem() << " " << tamanho_lista << endl;
	for (unsigned int i = 0; i < malha.r_num_elem(); i++)
	{
		Elemento elem = malha.r_elem(i);
		int num_nos = 0;
		if (elem.tipo == Tipo_Elem::Tri_Linear) num_nos = 3;
		if (elem.tipo == Tipo_Elem::Mini) num_nos = 3;
		if (elem.tipo == Tipo_Elem::Tri_Quadr) num_nos = 6;

		arquivo_saida << num_nos;
		for (int i = 0; i < num_nos; i++) arquivo_saida << " " << elem.nos[i];
		arquivo_saida << endl;
	}
	arquivo_saida << endl;

	arquivo_saida << "CELL_TYPES " << malha.r_num_elem() << endl;
	for (unsigned int i = 0; i < malha.r_num_elem(); i++)
	{
		Elemento el = malha.r_elem(i);
		if (el.tipo == Tipo_Elem::Tri_Linear) arquivo_saida << "5 ";
		if (el.tipo == Tipo_Elem::Mini) arquivo_saida << "5 ";
		if (el.tipo == Tipo_Elem::Tri_Quadr) arquivo_saida << "22 ";
	}
	arquivo_saida << endl << endl;

	arquivo_saida << "POINT_DATA " << vert << endl;
	arquivo_saida << "SCALARS vx double" << endl;
	arquivo_saida << "LOOKUP_TABLE default";
	for (unsigned int i = 0; i < vert; i++) arquivo_saida << endl << sol(i);
	arquivo_saida << endl << endl;

	arquivo_saida << "SCALARS vy double" << endl;
	arquivo_saida << "LOOKUP_TABLE default";
	for (unsigned int i = gl; i < gl + vert; i++) arquivo_saida << endl << sol(i);
	arquivo_saida << endl << endl;

	arquivo_saida << "SCALARS p double" << endl;
	arquivo_saida << "LOOKUP_TABLE default";
	for (unsigned int i = gl *2 ; i < gl * 2 + vert; i++) arquivo_saida << endl << sol(i);
	arquivo_saida << endl << endl;

	arquivo_saida.close();
}

//------ Método da Projeção ------

void NavierStokes2D::resolver_navier_stokes_projecao(double delta_t, unsigned long max_iter, double Re, string tag)
{
	iteracoes_por_saida = max_iter / max_saida;
	if (iteracoes_por_saida < 1) iteracoes_por_saida = 1;
	
	system_clock::time_point total1 = system_clock::now();
	
	cout << "Montando matrizes fixas K, M, Gx e Gy" << endl;
	system_clock::time_point t1 = system_clock::now();
	montar_matrizes_fixas(K, M, Gx, Gy);
	system_clock::time_point t2 = system_clock::now();
	cout << "A montagem das matrizes fixas levou " << duration_cast<seconds>(t2 - t1).count() << " segundos" << endl << endl;

	cout << "Montando B_proj, B_lumped_inv, G e D a partir das matrizes fixas" << endl;
	t1 = system_clock::now();	
	//B
	SparseMatrix<double, Eigen::RowMajor> B_proj(gl * 2, gl * 2);
	SparseMatrix<double> B_proj_col(gl * 2, gl * 2);
	vector<Triplet<double>> B_proj_triplets = montar_matriz_B_proj(true, delta_t, Re);
	B_proj.setFromTriplets(B_proj_triplets.begin(), B_proj_triplets.end());
	B_proj_col.setFromTriplets(B_proj_triplets.begin(), B_proj_triplets.end());
	
	//B_lumped
	SparseMatrix<double, Eigen::RowMajor> B_lumped_inversa(gl * 2, gl * 2);
	vector<Triplet<double>> B_lumped_triplets = montar_matriz_B_lumped_inversa(B_proj);
	B_lumped_inversa.setFromTriplets(B_lumped_triplets.begin(), B_lumped_triplets.end());
	//aplicar_cc_B_proj(B_lumped_inversa);

	//G
	SparseMatrix<double, Eigen::RowMajor> G(gl * 2, vert);
	SparseMatrix<double, Eigen::RowMajor> D(vert, gl * 2);
	SparseMatrix<double, Eigen::RowMajor> E(vert, vert);
	SparseMatrix<double, Eigen::RowMajor> D_B_inv_G(vert, vert);
	montar_matriz_G_e_D_proj(G, D);	
	
	t2 = system_clock::now();
	cout << "A montagem de B_proj, B_lumped_inv, G e D levou " << duration_cast<seconds>(t2 - t1).count() << " segundos" << endl << endl;

	cout << "Aplicando ccs nas matrizes B_proj, G e D" << endl;
	t1 = system_clock::now();
	aplicar_cc_B_proj(B_proj);
	aplicar_cc_D_e_G(D, G);
	aplicar_cc_E(E);
	//aplicar_cc_D_B_inv_G(D_B_inv_G);
	t2 = system_clock::now();

	D_B_inv_G = (D * B_lumped_inversa * G);
	D_B_inv_G = E + D_B_inv_G;

	cout << "A aplicação nas matrizes B_proj, G e D levou " << duration_cast<seconds>(t2 - t1).count() << " segundos" << endl << endl;

	cout << "Iniciando o semi Lagrangiano" << endl;
	t1 = system_clock::now();
	Semi_Lagrangiano SL(malha);
	t2 = system_clock::now();
	cout << "A inicializacao do semi Lagrangiano levou " << duration_cast<seconds>(t2 - t1).count() << " segundos" << endl << endl;

	cout << "Iniciando processo iterativo" << endl << endl;
	VectorXd vx(gl);
	VectorXd vy(gl);
	//VectorXd vx_novo(gl);
	//VectorXd vy_novo(gl);
	VectorXd v(gl * 2);
	VectorXd v_til(gl * 2);
	VectorXd p(vert);
	VectorXd M_v(gl * 2);
	VectorXd M_vx(gl);
	VectorXd M_vy(gl);
	VectorXd D_v_til(vert);
	VectorXd sol(gl * 2 + vert);
	v.fill(0);
	sol.fill(0);
	p.fill(0);

	for (unsigned int i = 0; i < max_iter; i++)
	{
		cout << "Iteracao: " << i << endl;
		t1 = system_clock::now();

		vx = v.head(gl);
		vy = v.segment(gl, gl);
		cout << "Utilizando semi lagrangiano" << endl;
		SL.semi_lagrangiano(malha, vx, vy, delta_t);
		SL.multiplicar_SL(vx);
		SL.multiplicar_SL(vy);		

		M_vx = M * ((double)1 / delta_t) * vx;
		M_vy = M * ((double)1 / delta_t) * vy;
		M_v.head(gl) = M_vx;
		M_v.segment(gl, gl) = M_vy;

		cout << "Aplicando condicoes de contorno ao vetor M_v" << endl;
		aplicar_cc_M_v(M_v);

		cout << "Resolvendo sistema linear B_proj * v_til = M_v" << endl << endl;
		Eigen::SparseLU<SparseMatrix<double>> solver;
		solver.compute(B_proj);
		if (solver.info() != EXIT_SUCCESS) cout << "Decomposicao falhou" << endl;
		v_til = solver.solve(M_v);

		D_v_til = D * v_til;
		aplicar_cc_D_v_til(D_v_til);

		cout << "Resolvendo sistema linear D_B_inv_G * p = D_v_til" << endl << endl;
		Eigen::SparseLU<SparseMatrix<double, Eigen::RowMajor>> solver2;
		solver2.compute(D_B_inv_G);
		if (solver2.info() != EXIT_SUCCESS) cout << "Decomposicao falhou" << endl;
		p = solver2.solve(D_v_til);

		v = v_til - (B_lumped_inversa * G * p);
		
		//sol só é usado para gerar a saída de dados, para não ter q criar outra funcao
		sol.head(gl*2) = v;
		sol.segment(gl * 2, vert) = p;
		gerar_arquivo_saida(sol, i, delta_t, Re, tag);
		t2 = system_clock::now();
		cout << "A iteracao " << i << " levou " << duration_cast<seconds>(t2 - t1).count() << " segundos" << endl << endl;
	}

	system_clock::time_point total2 = system_clock::now();
	cout << "A malha possui " << gl << " nos, e " << malha.r_num_elem() << " elementos." << endl;
	cout << "Todo procedimento (excluindo a leitura da malha) levou " << duration_cast<seconds>(total2 - total1).count() << " segundos" << endl;
}

void NavierStokes2D::aplicar_cc_B_proj(SparseMatrix<double, Eigen::RowMajor>& B_proj)
{
	unordered_map<unsigned long, double> cc_vx_map = malha.r_cc_vx();
	for (unordered_map<unsigned long, double>::iterator it_map = cc_vx_map.begin(); it_map != cc_vx_map.end(); it_map++)
	{
		for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(B_proj, it_map->first); it; ++it)
		{
			B_proj.coeffRef(it.row(), it.col()) = 0;
		}
		B_proj.coeffRef(it_map->first, it_map->first) = 1;
	}
	unordered_map<unsigned long, double> cc_vy_map = malha.r_cc_vy();
	for (unordered_map<unsigned long, double>::iterator it_map = cc_vy_map.begin(); it_map != cc_vy_map.end(); it_map++)
	{
		for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(B_proj, it_map->first + gl); it; ++it)
		{
			B_proj.coeffRef(it.row(), it.col()) = 0;
		}
		B_proj.coeffRef(it_map->first + gl, it_map->first + gl) = 1;
	}
}

void NavierStokes2D::aplicar_cc_D_e_G(SparseMatrix<double, Eigen::RowMajor>& D, SparseMatrix<double, Eigen::RowMajor>& G)
{
	unordered_map<unsigned long, double> cc_vx_map = malha.r_cc_vx();
	for (unordered_map<unsigned long, double>::iterator it_map = cc_vx_map.begin(); it_map != cc_vx_map.end(); it_map++)
	{
		for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(G, it_map->first); it; ++it)
		{
			G.coeffRef(it.row(), it.col()) = 0;
		}
		//G.coeffRef(it->first, it->first) = 1;
	}
	unordered_map<unsigned long, double> cc_vy_map = malha.r_cc_vy();
	for (unordered_map<unsigned long, double>::iterator it_map = cc_vy_map.begin(); it_map != cc_vy_map.end(); it_map++)
	{
		for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(G, it_map->first + gl); it; ++it)
		{
			G.coeffRef(it.row(), it.col()) = 0;
		}
		//G.coeffRef(it->first + gl, it->first) = 1;
	}
	unordered_map<unsigned long, double> cc_p_map = malha.r_cc_p();
	for (unordered_map<unsigned long, double>::iterator it_map = cc_p_map.begin(); it_map != cc_p_map.end(); it_map++)
	{
		for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(D, it_map->first); it; ++it)
		{
			D.coeffRef(it.row(), it.col()) = 0;
		}
		//D.coeffRef(it->first, it->first) = 1;
	}
}

void NavierStokes2D::aplicar_cc_E(SparseMatrix<double, Eigen::RowMajor>& E)
{
	unordered_map<unsigned long, double> cc_p_map = malha.r_cc_p();
	for (unordered_map<unsigned long, double>::iterator it_map = cc_p_map.begin(); it_map != cc_p_map.end(); it_map++)
	{
		E.coeffRef(it_map->first, it_map->first) = 1;
	}
}

void NavierStokes2D::aplicar_cc_M_v(VectorXd& M_v)
{
	unordered_map<unsigned long, double> cc_vx_map = malha.r_cc_vx();
	for (unordered_map<unsigned long, double>::iterator it_map = cc_vx_map.begin(); it_map != cc_vx_map.end(); it_map++)
	{
		M_v(it_map->first) = it_map->second;
	}
	unordered_map<unsigned long, double> cc_vy_map = malha.r_cc_vy();
	for (unordered_map<unsigned long, double>::iterator it_map = cc_vy_map.begin(); it_map != cc_vy_map.end(); it_map++)
	{
		M_v(it_map->first + gl) = it_map->second;
	}
}

void NavierStokes2D::aplicar_cc_D_v_til(VectorXd& D_v_til)
{
	unordered_map<unsigned long, double> cc_p_map = malha.r_cc_p();
	for (unordered_map<unsigned long, double>::iterator it_map = cc_p_map.begin(); it_map != cc_p_map.end(); it_map++)
	{
		D_v_til(it_map->first) = it_map->second;
	}
}

void NavierStokes2D::aplicar_cc_D_B_inv_G(SparseMatrix<double, Eigen::RowMajor>& D_B_inv_G)
{
	unordered_map<unsigned long, double> cc_p_map = malha.r_cc_p();
	for (unordered_map<unsigned long, double>::iterator it_map = cc_p_map.begin(); it_map != cc_p_map.end(); it_map++)
	{
		for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(D_B_inv_G, it_map->first); it; ++it)
		{
			D_B_inv_G.coeffRef(it.row(), it.col()) = 0;
		}
		//D_B_inv_G.coeffRef(it->first, it->first) = 1;
	}
}

vector<Triplet<double>> NavierStokes2D::montar_matriz_B_lumped_inversa(SparseMatrix<double, Eigen::RowMajor>& B)
{
	vector<Triplet<double>> triplets;
	for (int i = 0; i < B.outerSize(); ++i)
	{
		double valor = 0;
		for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(B, i); it; ++it)
		{
			valor += it.value();
		}
		triplets.push_back(Triplet<double>(i, i, (double)1/valor));
	}
	return triplets;
}

vector<Triplet<double>> NavierStokes2D::montar_matriz_B_proj(bool transiente, double delta_t, double Re)
{
	vector<Triplet<double>> triplets;
	for (int i = 0; i < K.outerSize(); ++i)
	{
		for (SparseMatrix<double>::InnerIterator it(K, i); it; ++it)
		{
			triplets.push_back(Triplet<double>((long)it.row(), (long)it.col(), it.value() / Re));
			triplets.push_back(Triplet<double>((long)it.row() + gl, (long)it.col() + gl, it.value() / Re));
		}
	}
	if (transiente == true)
	{
		for (int i = 0; i < M.outerSize(); ++i)
		{
			for (SparseMatrix<double>::InnerIterator it(M, i); it; ++it)
			{
				triplets.push_back(Triplet<double>((long)it.row(), (long)it.col(), it.value() / delta_t));
				triplets.push_back(Triplet<double>((long)it.row() + gl, (long)it.col() + gl, it.value() / delta_t));
			}
		}
	}
	return triplets;
}

void NavierStokes2D::montar_matriz_G_e_D_proj(SparseMatrix<double, Eigen::RowMajor>& G, SparseMatrix<double, Eigen::RowMajor>& D)
{
	vector<Triplet<double>> G_triplets;
	vector<Triplet<double>> D_triplets;
	for (int i = 0; i < Gx.outerSize(); ++i)
	{
		for (SparseMatrix<double>::InnerIterator it(Gx, i); it; ++it)
		{
			G_triplets.push_back(Triplet<double>((long)it.row(), (long)it.col(), it.value()));
			D_triplets.push_back(Triplet<double>((long)it.col(), (long)it.row(), it.value()));
		}
	}
	for (int i = 0; i < Gy.outerSize(); ++i)
	{
		for (SparseMatrix<double>::InnerIterator it(Gy, i); it; ++it)
		{
			G_triplets.push_back(Triplet<double>((long)it.row() + gl, (long)it.col(), it.value()));
			D_triplets.push_back(Triplet<double>((long)it.col(), (long)it.row() + gl, it.value()));
		}
	}

	G.setFromTriplets(G_triplets.begin(), G_triplets.end());
	D.setFromTriplets(D_triplets.begin(), D_triplets.end());
}

//------ CC na linha e na coluna ------

//void NavierStokes2D::aplicar_cc_linha_e_coluna(SparseMatrix<double, Eigen::RowMajor>& A)
//{
//	unordered_map<unsigned long, CC> cc_vx_map = malha.r_cc_vx();
//	for (unordered_map<unsigned long, CC>::iterator it = cc_vx_map.begin(); it != cc_vx_map.end(); it++)
//	{
//		for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(A, it->first); it; ++it)
//		{
//			A.coeffRef(it.row(), it.col()) = 0;
//		}
//		A.coeffRef(it->first, it->first) = 1;
//	}
//	unordered_map<unsigned long, CC> cc_vy_map = malha.r_cc_vy();
//	for (unordered_map<unsigned long, CC>::iterator it = cc_vy_map.begin(); it != cc_vy_map.end(); it++)
//	{
//		for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(A, it->first + gl); it; ++it)
//		{
//			A.coeffRef(it.row(), it.col()) = 0;
//		}
//		A.coeffRef(it->first + gl, it->first + gl) = 1;
//	}
//	unordered_map<unsigned long, CC> cc_p_map = malha.r_cc_p();
//	for (unordered_map<unsigned long, CC>::iterator it = cc_p_map.begin(); it != cc_p_map.end(); it++)
//	{
//		for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(A, it->first + gl * 2); it; ++it)
//		{
//			A.coeffRef(it.row(), it.col()) = 0;
//		}
//		A.coeffRef(it->first + gl * 2, it->first + gl * 2) = 1;
//	}
//
//	for (int k = 0; k < A.outerSize(); ++k) //Percorre todas as linhas. De 0 a k
//	{
//		vector<unsigned long> indices_col;
//		
//		//it(A, k) -> Percorra a linha k da matriz A
//		for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(A, k); it; ++it)
//		{
//			//cout << it.row() << " " << it.col() << endl;
//			indices_col.push_back(it.col());
//		}
//
//		for (unsigned long i = 0; i < indices_col.size(); i++)
//		{			
//			unordered_map<unsigned long, CC>::iterator it_vx_map = cc_vx_map.find(indices_col[i]);
//			if (it_vx_map != cc_vx_map.end())
//			{
//				A.coeffRef(k, indices_col[i]) = 0;
//			}
//			unordered_map<unsigned long, CC>::iterator it_vy_map = cc_vy_map.find(indices_col[i] - gl);
//			if (it_vy_map != cc_vy_map.end())
//			{
//				A.coeffRef(k, indices_col[i]) = 0;
//			}
//			unordered_map<unsigned long, CC>::iterator it_p_map = cc_p_map.find(indices_col[i] - gl * 2);
//			if (it_p_map != cc_p_map.end())
//			{
//				A.coeffRef(k, indices_col[i]) = 0;
//			}
//		}
//	}
//}

void NavierStokes2D::gerar_linha_prop(int num_pontos, double xo, double yo, double x, double y)
{
	linha_prop_x = x;
	linha_prop_xo = xo;
	linha_prop_y = y;
	linha_prop_yo = yo;
	linha_prop_num_pontos = num_pontos;
}

void NavierStokes2D::avaliar_prop(VectorXd& prop, Semi_Lagrangiano& sl, string tag, int num_pontos, double xo, double yo, double x, double y)
{
	if (num_pontos == 0) return;
	
	fstream saida_vx;
	fstream saida_vy;
	fstream saida_p;
	saida_vx.open(malha.r_nome() + " " + tag + "_vx.txt", ios::out);
	saida_vy.open(malha.r_nome() + " " + tag + "_vy.txt", ios::out);
	saida_p.open(malha.r_nome() + " " + tag + "_p.txt", ios::out);

	double lx = (x - xo) / (num_pontos - 1);
	double ly = (y - yo) / (num_pontos - 1);

	for (unsigned int i = 0; i < num_pontos; i++)
	{
		double px = lx * i + xo;
		double py = ly * i + yo;

		unsigned el_index = sl.busca_linear(malha, px, py);
		if (el_index == -1) cout << "PROBLEMA AVALIAR PROP" << endl;

		Elemento elem = malha.r_elem(el_index);
		No no1 = malha.r_no(elem.nos[0]);
		No no2 = malha.r_no(elem.nos[1]);
		No no3 = malha.r_no(elem.nos[2]);
		double peso1 = ((no2.y - no3.y) * (px - no3.x) + (no3.x - no2.x) * (py - no3.y)) /
			((no2.y - no3.y) * (no1.x - no3.x) + (no3.x - no2.x) * (no1.y - no3.y));
		double peso2 = ((no3.y - no1.y) * (px - no3.x) + (no1.x - no3.x) * (py - no3.y)) /
			((no2.y - no3.y) * (no1.x - no3.x) + (no3.x - no2.x) * (no1.y - no3.y));
		double peso3 = 1 - peso1 - peso2;

		double vx_x = sqrt(pow(lx * i, 2) + pow(ly * i, 2));
		double vx_y =  prop[no1.id] * peso1 + prop[no2.id] * peso2 + prop[no3.id] * peso3;
		saida_vx << vx_x << " " << vx_y << endl;

		double vy_x = sqrt(pow(lx * i, 2) + pow(ly * i, 2));
		double vy_y = prop[no1.id + gl] * peso1 + prop[no2.id+gl] * peso2 + prop[no3.id+gl] * peso3;
		saida_vy << vy_x << " " << vy_y << endl;

		double p_x = sqrt(pow(lx * i, 2) + pow(ly * i, 2));
		double p_y = prop[no1.id + gl*2] * peso1 + prop[no2.id + gl*2] * peso2 + prop[no3.id + gl*2] * peso3;
		saida_p << p_x << " " << p_y << endl;
	}
}