#include "NavierStokes2D.h"
#include "NS_Bifasico.h"

int main()
{
	//NavierStokes2D n(".\\dreno\\dreno");
	//NavierStokes2D n(".\\Backward_facing_step\\backward_facing_step");
	//NavierStokes2D n(".\\Poiseuille\\poiseuille4");
	//NavierStokes2D n(".\\malha_teste\\malha_teste");
	//NavierStokes2D n(".\\Multi_furo\\multi_furo");
	//NavierStokes2D n(".\\Liddriven\\liddriven");
	NS_Bifasico nb(".\\Bifasico\\Reynolds1000_estatico\\caixa2_2", ".\\Bifasico\\Reynolds1000_estatico\\bolha_circular");
	
	double Reynolds = 1000;
	double dt = 0.001;
	long max_iter = 2000;

	double mi_in = 0.01;
	double rho_in = 0.001;
	double mi_out = 1;
	double rho_out = 1;
	
	//n.resolver_stokes_permanente();	
	//n.resolver_stokes_transiente(1, 30);
	//n.gerar_linha_prop(100, 1, 0, 1, 1);
	//n.resolver_navier_stokes(dt, max_iter, Reynolds, "oi");
	nb.gerar_linha_prop(400, 0, 1, 2, 1);
	nb.resolver_navier_stokes(dt, max_iter, Reynolds, mi_in, rho_in, mi_out, rho_out, "11_maio");
	//n.resolver_navier_stokes_projecao(dt, max_iter, Reynolds, "proj");
}