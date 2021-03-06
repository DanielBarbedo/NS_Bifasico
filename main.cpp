#include "NavierStokes2D.h"
#include "NS_Bifasico.h"

int main()
{
	//NavierStokes2D n(".\\Bifasico\\poiseuille_vert");
	//NavierStokes2D n(".\\Backward_facing_step\\backward_facing_step");
	//NavierStokes2D n(".\\Tubo_furo\\tubo_furo2");
	//NavierStokes2D n(".\\Cilindro2\\cilindro2");
	//NavierStokes2D n(".\\Malha_teste\\malha_teste");
	//NavierStokes2D n(".\\Liddriven\\liddriven");
	NS_Bifasico n(".\\Bifasico\\poiseuille_vert", ".\\Bifasico\\bolha");
	
	//n.resolver_stokes_permanente();	
	//n.resolver_stokes_transiente(1, 30);
	n.resolver_navier_stokes(0.01, 60, 0.01, 0.001, 1, 1, "bi");
	//n.resolver_navier_stokes_projecao(0.01, 30, 1, 1);

	//Malha m;
	//m.ler_malha(".\\Tubo_furo\\tubo_furo");

	//VectorXd vx(m.r_num_nos());
	//VectorXd vy(m.r_num_nos());
	//vx.fill(1);
	//vy.fill(0);

	//Semi_Lagrangiano sl(m);
	//sl.semi_lagrangiano(m, vx, vy, 0.2, vx, false);

	//Malha_Bolha m;
	//m.ler_malha(".\\Bifasico\\bolha");
	//m.mostrar_malha();
}