#include "NS_Bifasico.h"

NS_Bifasico::NS_Bifasico(string nome, string nome_bolha)
{
	malha.ler_malha(nome);
	malha.gerar_elemento_mini();
	//malha.mostrar_malha();
	malha_bolha.ler_malha(nome_bolha);

	gl = malha.r_num_nos();
	vert = malha.r_num_vertices();

	cout << "--- Malha lida, com " << gl << " nos e " << malha.r_num_elem() << " elementos. ---" << endl << endl;
}

void NS_Bifasico::resolver_navier_stokes(double delta_t, unsigned long max_iter, double mi_in, double rho_in, double mi_out, double rho_out, string tag)
{	
	cout << "Alocando matrizes e vetores" << endl;
	SparseMatrix<double> K(gl, gl);
	SparseMatrix<double> M(gl, gl);
	SparseMatrix<double> Gx(gl, vert);
	SparseMatrix<double> Gy(gl, vert);
	SparseMatrix<double, Eigen::RowMajor> A(2 * gl + vert, 2 * gl + vert);
	VectorXd B(gl * 2 + vert);
	VectorXd sol(gl * 2 + vert);
	VectorXd vx(gl);
	VectorXd vy(gl);
	VectorXd vx_temp(gl);
	VectorXd vy_temp(gl);
	VectorXd rho_vec(vert);
	vx.fill(0);
	vy.fill(0);
	sol.fill(0);

	cout << "Iniciando Semi Lagrangiano" << endl;
	Semi_Lagrangiano sl(malha);

	for (unsigned int i = 0; i < max_iter; i++)
	{
		cout << endl << "ITERACAO: " << i << endl << endl;
		
		cout << "Montando matrizes" << endl;
		montar_matrizes(K, M, Gx, Gy, mi_in, rho_in, mi_out, rho_out, sl, rho_vec);
		montar_matriz_A(K, M, Gx, Gy, A, delta_t, 1, 1); //------- MUDAR AQUI DEPOIS, NI E RHO PRA REYNOLDS ------

		cout << "Atualizando posicao da bolha" << endl;
		malha_bolha.atualizar_posicao(vx, vy, delta_t);		

		vx_temp = sol.head(gl);
		vy_temp = sol.segment(gl, gl);
		sl.semi_lagrangiano(malha, vx, vy, delta_t, vx_temp, false);
		sl.semi_lagrangiano(malha, vx, vy, delta_t, vy_temp, false);
		vx = vx_temp;
		vy = vy_temp;

		B.fill(0);
		B.head(gl) = M * ((double)1 / delta_t) * vx;
		B.segment(gl, gl) = M * ((double)1 / delta_t) * vy;

		cout << "Aplicando ccs" << endl;
		aplicar_cc_A(A);
		aplicar_cc_B(B);

		cout << "Resolvendo sistema linear" << endl << endl;
		//Eigen::BiCGSTAB<SparseMatrix<double>> solver;
		Eigen::SparseLU<SparseMatrix<double>> solver;
		solver.compute(A);
		if (solver.info() != EXIT_SUCCESS) cout << "Decomposicao falhou" << endl;
		sol = solver.solve(B);

		gerar_arquivo_saida(sol, rho_vec, i, delta_t, 1, 1, tag); //------- MUDAR AQUI DEPOIS, NI E RHO PRA REYNOLDS ------
		gerar_saida_bolha(i, tag);
	}
}

void NS_Bifasico::gerar_saida_bolha(unsigned long iter, string tag)
{
	fstream arquivo_saida;
	//arquivo_saida.open(malha.r_nome() + "_ni_" + to_string(ni) + "_rho_" + to_string(rho) 
	//	+ "_dt_" + to_string(delta_t) + "_" + tag + to_string(iter) + ".vtk", ios::out);
	arquivo_saida.open(malha.r_nome() + "_bolha_" + "_" + tag + "_" + to_string(iter) + ".vtk", ios::out);
	unsigned long gl = malha_bolha.r_num_nos();

	arquivo_saida << "# vtk DataFile Version 1.0" << endl;
	arquivo_saida << "Funcao Navier-Stokes 2D" << endl;
	arquivo_saida << "ASCII" << endl;
	arquivo_saida << "DATASET UNSTRUCTURED_GRID" << endl;
	arquivo_saida << endl;

	arquivo_saida << "POINTS " << gl << " double" << endl;
	for (unsigned int i = 0; i < gl; i++)
	{
		No no = malha_bolha.r_no(i);
		arquivo_saida << no.x << " " << no.y << " " << no.z << endl;
	}
	arquivo_saida << endl;

	unsigned long tamanho_lista = malha_bolha.r_num_nos() * 2;
	arquivo_saida << "CELLS " << malha_bolha.r_num_nos() << " " << tamanho_lista << endl;
	for (unsigned int i = 0; i < malha_bolha.r_num_nos(); i++)
	{
		arquivo_saida << 1 << " " << i << endl;
	}
	arquivo_saida << endl;

	arquivo_saida << "CELL_TYPES " << malha_bolha.r_num_nos() << endl;
	for (unsigned int i = 0; i < malha_bolha.r_num_nos(); i++) arquivo_saida << "1 ";
	arquivo_saida << endl << endl;
}

void NS_Bifasico::gerar_arquivo_saida(const VectorXd& sol, const VectorXd& rho_vec, unsigned long iter, double delta_t, double ni, double rho, string tag)
{
	fstream arquivo_saida;
	//arquivo_saida.open(malha.r_nome() + "_ni_" + to_string(ni) + "_rho_" + to_string(rho) 
	//	+ "_dt_" + to_string(delta_t) + "_" + tag + to_string(iter) + ".vtk", ios::out);
	arquivo_saida.open(malha.r_nome() + "_dt_" + to_string(delta_t) + "_ni_" + to_string(ni) + "_" + tag + "_"
		+ to_string(iter) + ".vtk", ios::out);
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
	for (unsigned int i = gl * 2; i < gl * 2 + vert; i++) arquivo_saida << endl << sol(i);
	arquivo_saida << endl << endl;

	arquivo_saida << "SCALARS rho double" << endl;
	arquivo_saida << "LOOKUP_TABLE default";
	for (unsigned int i = 0; i < vert; i++) arquivo_saida << endl << rho_vec(i);
	arquivo_saida << endl << endl;

	arquivo_saida.close();
}

void NS_Bifasico::aplicar_cc_A(SparseMatrix<double, Eigen::RowMajor>& A)
{
	unordered_map<unsigned long, double> cc_vx_map = malha.r_cc_vx();
	for (unordered_map<unsigned long, double>::iterator it = cc_vx_map.begin(); it != cc_vx_map.end(); it++)
	{
		for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(A, it->first); it; ++it)
		{
			//it.valueRef() = 0;
			A.coeffRef(it.row(), it.col()) = 0;
		}
		A.coeffRef(it->first, it->first) = 1;
	}
	unordered_map<unsigned long, double> cc_vy_map = malha.r_cc_vy();
	for (unordered_map<unsigned long, double>::iterator it = cc_vy_map.begin(); it != cc_vy_map.end(); it++)
	{
		for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(A, it->first + gl); it; ++it)
		{
			//it.valueRef() = 0;
			A.coeffRef(it.row(), it.col()) = 0;
		}
		A.coeffRef(it->first + gl, it->first + gl) = 1;
	}
	unordered_map<unsigned long, double> cc_p_map = malha.r_cc_p();
	for (unordered_map<unsigned long, double>::iterator it = cc_p_map.begin(); it != cc_p_map.end(); it++)
	{
		for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(A, it->first + gl * 2); it; ++it)
		{
			//it.valueRef() = 0;
			A.coeffRef(it.row(), it.col()) = 0;
		}
		A.coeffRef(it->first + gl * 2, it->first + gl * 2) = 1;
	}
}

void NS_Bifasico::aplicar_cc_B(VectorXd& B)
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

void NS_Bifasico::montar_matriz_A(SparseMatrix<double>& K, SparseMatrix<double>& M, SparseMatrix<double>& Gx, SparseMatrix<double>& Gy, 
	SparseMatrix<double, Eigen::RowMajor>& A, double delta_t, double ni, double rho)
{
	vector<Triplet<double>> A_triplets;
	for (int i = 0; i < K.outerSize(); ++i)
	{
		for (SparseMatrix<double>::InnerIterator it(K, i); it; ++it)
		{
			//it.index(); // inner index, here it is equal to it.row()
			A_triplets.push_back(Triplet<double>((long)it.row(), (long)it.col(), ni * it.value()));
			A_triplets.push_back(Triplet<double>((long)it.row() + gl, (long)it.col() + gl, ni * it.value()));
		}
	}
	for (int i = 0; i < M.outerSize(); ++i)
	{
		for (SparseMatrix<double>::InnerIterator it(M, i); it; ++it)
		{
			A_triplets.push_back(Triplet<double>((long)it.row(), (long)it.col(), it.value() / delta_t));
			A_triplets.push_back(Triplet<double>((long)it.row() + gl, (long)it.col() + gl, it.value() / delta_t));
		}
	}
	for (int i = 0; i < Gx.outerSize(); ++i)
	{
		for (SparseMatrix<double>::InnerIterator it(Gx, i); it; ++it)
		{
			A_triplets.push_back(Triplet<double>((long)it.row(), (long)it.col() + gl * 2, ((double)1 / rho) * it.value()));
			A_triplets.push_back(Triplet<double>((long)it.col() + gl * 2, (long)it.row(), -it.value()));
		}
	}
	for (int i = 0; i < Gy.outerSize(); ++i)
	{
		for (SparseMatrix<double>::InnerIterator it(Gy, i); it; ++it)
		{
			A_triplets.push_back(Triplet<double>((long)it.row() + gl, (long)it.col() + gl * 2, ((double)1 / rho) * it.value()));
			A_triplets.push_back(Triplet<double>((long)it.col() + gl * 2, (long)it.row() + gl, -it.value()));
		}
	}

	A.setFromTriplets(A_triplets.begin(), A_triplets.end());
}

void NS_Bifasico::montar_matrizes(SparseMatrix<double>& K, SparseMatrix<double>& M, SparseMatrix<double>& Gx, SparseMatrix<double>& Gy,
	double mi_in, double rho_in, double mi_out, double rho_out, Semi_Lagrangiano& sl, VectorXd& rho_vec)
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
		No no4 = malha.r_no(elem.nos[3]); //No do centro
		
		double rho_no1 = rho_in * calc_heavyside(sl, no1.id) + rho_out * (1 - calc_heavyside(sl, no1.id));
		double rho_no2 = rho_in * calc_heavyside(sl, no2.id) + rho_out * (1 - calc_heavyside(sl, no2.id));
		double rho_no3 = rho_in * calc_heavyside(sl, no3.id) + rho_out * (1 - calc_heavyside(sl, no3.id));
		double rho = (rho_no1 + rho_no2 + rho_no3) / 3;

		//cout << rho << endl;

		rho_vec[no1.id] = rho_no1;
		rho_vec[no2.id] = rho_no2;
		rho_vec[no3.id] = rho_no3;
		
		double mi_no1 = mi_in * calc_heavyside(sl, no1.id) + mi_out * (1 - calc_heavyside(sl, no1.id));
		double mi_no2 = mi_in * calc_heavyside(sl, no2.id) + mi_out * (1 - calc_heavyside(sl, no2.id));
		double mi_no3 = mi_in * calc_heavyside(sl, no3.id) + mi_out * (1 - calc_heavyside(sl, no3.id));
		double mi = (mi_no1 + mi_no2 + mi_no3) / 2;

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
		//------ ADICAO BIFASICO ------
		mat_ke_mini = mat_ke_mini * mi;

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
		//------ ADICAO BIFASICO ------
		mat_me_mini = mat_me_mini * rho;

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

double NS_Bifasico::calc_heavyside(Semi_Lagrangiano& sl, unsigned long no_ref_index)
{
	No no_ref = malha.r_no(no_ref_index);	
	Elemento el = malha_bolha.r_elem(0);

	No no1 = malha_bolha.r_no(el.nos[0]);
	No no2 = malha_bolha.r_no(el.nos[1]);
	double x_el = (no1.x + no2.x) / 2;
	double y_el = (no1.y + no2.y) / 2;

	//No no_ref;
	//no_ref.id = 10;
	//no_ref.x = 0.3;
	//no_ref.y = 0.55;
	//no_ref.z = 0;

	double dist = pow(no_ref.x - x_el, 2) + pow(no_ref.y - y_el, 2); //Como é comparação, dispensa-se a raiz
	long index_menor_dist = 0;

	for (unsigned int i = 1; i < malha_bolha.r_num_nos(); i++)
	{
		Elemento el = malha_bolha.r_elem(i);
		No no1 = malha_bolha.r_no(el.nos[0]);
		No no2 = malha_bolha.r_no(el.nos[1]);
		double x_el = (no1.x + no2.x) / 2;
		double y_el = (no1.y + no2.y) / 2;

		double novo_dist = pow(no_ref.x - x_el, 2) + pow(no_ref.y - y_el, 2); //Como é comparação, dispensa-se a raiz
		if (novo_dist < dist)
		{
			dist = novo_dist;
			index_menor_dist = i;
		}
	}

	el = malha_bolha.r_elem(index_menor_dist);
	no1 = malha_bolha.r_no(el.nos[0]);
	no2 = malha_bolha.r_no(el.nos[1]);
	
	double px = (no2.x + no1.x) / 2;
	double py = (no2.y + no1.y) / 2;
	//Você pode adicionar depois uma checagem pro sentido de giro do círculo.
	//A depender da sequência de pontos (horário ou anti-horário) os vetores dos
	//elementos vão apontar pra fora ou pra dentro.

	double v_ref_x = no_ref.x - px;
	double v_ref_y = no_ref.y - py;
	double dx = no2.x - no1.x;
	double dy = no2.y - no1.y;
	double vetor_x = dy;
	double vetor_y = -dx;

	double escalar = v_ref_x * vetor_x + v_ref_y * vetor_y;
	if (escalar < 0) return 1;
	else return 0;
}