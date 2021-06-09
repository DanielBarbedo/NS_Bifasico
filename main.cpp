#include "NavierStokes2D.h"
#include "NS_Bifasico.h"

int main()
{
	//NavierStokes2D n(".\\dreno\\dreno");
	//NavierStokes2D n(".\\Backward_facing_step\\backward_facing_step");
	//NavierStokes2D n(".\\Poiseuille\\poiseuille3");
	//NavierStokes2D n(".\\malha_teste\\malha_teste");
	//NavierStokes2D n(".\\Multi_furo\\multi_furo");
	//NavierStokes2D n(".\\Liddriven\\liddriven");
	//NS_Bifasico nb(".\\Bi\\Reynolds1000_5k_19_05\\caixa2_2", ".\\Bi\\Reynolds1000_5k_19_05\\bolha_elipse");
	NS_Bifasico nb(".\\Bi\\Reynolds100\\caixa2_2", ".\\Bi\\Reynolds100\\bolha_circular");
	
	double Reynolds = 100;
	double dt = 0.001;
	long max_iter = 2000;

	//bolha
	//double mi_in = 0.01;
	//double rho_in = 0.001;
	//double mi_out = 1;
	//double rho_out = 1;
	
	//gota
	//double mi_in = 1;
	//double rho_in = 1;
	//double mi_out = 0.01;
	//double rho_out = 0.001;

	double mi_in = 0.01;
	double rho_in = 0.001;
	double mi_out = 1;
	double rho_out = 1;
	double Weber = 20;

	double l_interface = 0.2;
	
	//n.resolver_stokes_permanente();	
	//n.resolver_stokes_transiente(1, 30);
	//n.gerar_linha_prop(100, 1, 0, 1, 1);
	//n.resolver_navier_stokes(dt, max_iter, Reynolds, "oi");
	nb.alterar_numero_arquivos_saida(2000);
	nb.gerar_linha_prop(400, 0, 1, 2, 1);
	nb.resolver_navier_stokes(dt, max_iter, Reynolds, Weber, mi_in, rho_in, mi_out, rho_out, l_interface, "gota_circ_grav");
	//n.resolver_navier_stokes_projecao(dt, max_iter, Reynolds, "proj");
}