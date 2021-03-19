#include "NavierStokes2D.h"
#include "NS_Bifasico.h"

int main()
{
	//NavierStokes2D n(".\\Bifasico\\poiseuille_vert");
	//NavierStokes2D n(".\\Backward_facing_step\\backward_facing_step");
	//NavierStokes2D n(".\\Tubo_furo\\tubo_furo2");
	//NavierStokes2D n(".\\Cilindro2\\cilindro2");
	NavierStokes2D n(".\\Rebaixo\\rebaixo");
	//NavierStokes2D n(".\\Liddriven\\liddriven");
	//NS_Bifasico n(".\\Bifasico\\poiseuille_vert_no_inflow", ".\\Bifasico\\bolha_elipse");
	
	//n.resolver_stokes_permanente();	
	//n.resolver_stokes_transiente(1, 30);
	n.resolver_navier_stokes(0.01, 100, 0.001, 1, false, "ns");
	//n.resolver_navier_stokes(0.01, 100, 1, 1, 0.01, 0.001, "bi");
	//n.resolver_navier_stokes_projecao(0.01, 30, 1, 1);
}