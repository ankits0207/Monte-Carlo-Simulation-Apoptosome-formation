package com.sssbb.assignment1;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Random;

public class Monte_Carlo_Simulation_cytoc_apaf {
	double delta_t = Math.pow(10, -4);
	// static double mc_steps = Math.pow(10, 8);
	static double mc_steps = Math.pow(10, 6);
	double total_time = delta_t*mc_steps;
	static int cytoc = 100;
	static int apaf = 50;
	static int lattice_dim = 60;
	static int total_molecules = cytoc + apaf;
	static double[] p_rwalk = {1/6,2/6,3/6,4/6,5/6,1};
	static double p_apoptosome_form = 2.8*Math.pow(10, -5);
	static double p_apoptosome_break = 5.7*Math.pow(10, -7);
	double p_dimer_form = 1;
	double p_dimer_break = 0;
	static double p_apoptosome_diffusion = 0.1;
	static int apoptosome_count = 0;
	static int dimer_count = 0;
	static int[] global_apoptosome_count = new int[(int) mc_steps];
	static int[] global_dimer_count = new int[(int) mc_steps];
	static int my_object_array_index = 0;
	
	static int flag = 0;
	static int flag2 = 0;
	
	static Molecule_Class[][][] lattice = new Molecule_Class[lattice_dim][lattice_dim][lattice_dim];
	static Molecule_Class[] my_object_array = new Molecule_Class[total_molecules];
		
	public static void initialization()
	{
		my_object_array_index = 0;
//		for (int i = 0; i < my_object_array.length; i++)
//		{
//			my_object_array[i] = null;
//		}
//		for (int i = 0; i < lattice.length; i++)
//		{
//			lattice[i] = null;
//		}
		//my_object_array = null;
		//lattice = null;
		for (int i = 0; i < (int) mc_steps; i++) {
			global_apoptosome_count[i] = 0;
			global_dimer_count[i] = 0;
		}
		for (int i = 0; i < cytoc; i++)
		{
			Random rand = new Random();
			int r1 = rand.nextInt(lattice_dim);
			int r2 = rand.nextInt(lattice_dim);
			int r3 = rand.nextInt(lattice_dim);
			while(lattice[r1][r2][r3] != null)
			{
				r1 = rand.nextInt(lattice_dim);
				r2 = rand.nextInt(lattice_dim);
				r3 = rand.nextInt(lattice_dim);
			}
			Molecule_Class mol_obj = new Molecule_Class(1, 0, r1, r2, r3, 0);
			lattice[r1][r2][r3] = mol_obj;
			my_object_array[my_object_array_index] = mol_obj;
			my_object_array_index ++;
		}
		for (int j = 0; j < apaf; j++)
		{
			Random rand = new Random();
			int r1 = rand.nextInt(lattice_dim);
			int r2 = rand.nextInt(lattice_dim);
			int r3 = rand.nextInt(lattice_dim);
			while(lattice[r1][r2][r3] != null)
			{
				r1 = rand.nextInt(lattice_dim);
				r2 = rand.nextInt(lattice_dim);
				r3 = rand.nextInt(lattice_dim);
			}
			Molecule_Class mol_obj = new Molecule_Class(2, 0, r1, r2, r3, 0);
			lattice[r1][r2][r3] = mol_obj;
			my_object_array[my_object_array_index] = mol_obj;
			my_object_array_index ++;
		}
		
	}
	
	public static void diffusion_cytoc_apaf(int r)
	{
		// By: Divyanshu Srivastava
		// This function randomly picks a molecule in the lattice, and tries to diffuse 
		// it in its surroundings.
		double rand = Math.random();
		int x, y, z;
		x = my_object_array[r].xaxis;
        y = my_object_array[r].yaxis;
        z = my_object_array[r].zaxis;
        
		if (my_object_array[r].molecule_bonded_with == 0) {
			// Only free molecules are allowed to move
            
            if (x != 0 && rand < p_rwalk[0] && lattice[x-1][y][z] == null) {
            	// move the molecule one step left        
	            // update lattice
	            lattice[x-1][y][z] = my_object_array[r];
                lattice[x][y][z] = null;
                // update my_object_array
                my_object_array[r].xaxis = x-1;
            }
            else if (x != (lattice_dim-1) && rand < p_rwalk[1] && lattice[x+1][y][z] == null) {
                // move the molecule one step right
                // update lattice
                lattice[x+1][y][z] = my_object_array[r];
            	lattice[x][y][z] = null;	
                // update my_object_array
                my_object_array[r].xaxis = x+1;
            }
            else if (y != 0 && rand < p_rwalk[2] && lattice[x][y-1][z] == null) {
                // move the molecule one step up
                // update lattice
                lattice[x][y-1][z] = my_object_array[r];
                lattice[x][y][z] = null;	
                // update my_object_array
                my_object_array[r].yaxis = y-1;
            }
            else if (y != (lattice_dim-1) && rand < p_rwalk[3] && lattice[x][y+1][z] == null) {
                // move the molecule one step down
                // update lattice
                lattice[x][y+1][z] = my_object_array[r];
                lattice[x][y][z] = null;
                // update my_object_array
                my_object_array[r].yaxis = y+1;
            }
            else if (z != 0 && rand < p_rwalk[4] && lattice[x][y][z-1] == null) {
                // move the molecule one step front
                // update lattice
                lattice[x][y][z-1] = my_object_array[r];
                lattice[x][y][z] = null;
                // update my_object_array
                my_object_array[r].zaxis = z-1;
            }
            else if (z != (lattice_dim-1) && rand < p_rwalk[5] && lattice[x][y][z+1] == null) {
                // move the molecule one step back
                // update lattice
                lattice[x][y][z+1] = my_object_array[r];
                lattice[x][y][z] = null;
                // update my_object_array
                my_object_array[r].zaxis = z+1;
            }
		}
	}
	
	public static void diffusion_apoptosome(int r) {
		// By: Divyanshu Srivastava
		
		if(my_object_array[r].molecule_type == 3) {
			// checking if the molecule is apoptosome
			Random rand = new Random();
			int choice = rand.nextInt(6);	// for six possible directions
			choice = choice+1;
			int x = my_object_array[r].xaxis;
			int y = my_object_array[r].yaxis;
			int z = my_object_array[r].zaxis;
			int partner_x = x;
			int partner_y = y;
			int partner_z = z;
			switch (my_object_array[r].bond_direction) {
			case (1):partner_x = x+1;
				break;
			case (2):partner_x = x-1;
				break;
			case (3):partner_y = y+1;
				break;
			case (4):partner_y = y-1;
				break;
			case (5):partner_z = z+1;
				break;
			case (6):partner_z = z-1;
				break;
			default:break;
			
			}
			
			int t;	// stores index of partner
			for (t = 0; t<total_molecules; t++) {
				if ((my_object_array[t].xaxis == partner_x) && (my_object_array[t].yaxis == partner_y) && (my_object_array[t].zaxis == partner_z)) {
//					System.out.println("Before BREAK" );
					break;
				}
			}
			if (t == 150) {
				System.out.println("Code broken @ t = 150");
				System.out.println("r = " + r);
				System.out.println("molecule location");
				System.out.println(x);
				System.out.println(y);
				System.out.println(z);
				System.out.println("partner location");
				System.out.println(partner_x);
				System.out.println(partner_y);
				System.out.println(partner_z);
				for (int lk = 0; lk<total_molecules; lk++) {
					System.out.println("Molecule: " + lk);
					System.out.println(my_object_array[lk].xaxis);
					System.out.println(my_object_array[lk].yaxis);
					System.out.println(my_object_array[lk].zaxis);
					System.out.println(my_object_array[lk].bond_direction);
				}
				try {
					System.in.read();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			if (choice == 1) {
				// Right
				if (choice != my_object_array[r].bond_direction  && x<(lattice_dim-1) && lattice[x+1][y][z] == null) {
					if (my_object_array[r].bond_direction == 2) {
						// Same bond orientation in opposite direction
						// All set to move.
						// lattice update
						lattice[x+1][y][z] = my_object_array[r];
						lattice[x][y][z] = lattice[partner_x][partner_y][partner_z];
						lattice[partner_x][partner_y][partner_z] = null;
						// my_object_array update
						my_object_array[r].xaxis = x+1;
						my_object_array[t].xaxis = partner_x+1;
					} else {
						// Different bond orientation
						if (lattice[partner_x+1][partner_y][partner_z] == null) {
							// All set to move
							// lattice update
							lattice[x+1][y][z] = my_object_array[r];
							lattice[x][y][z] = null;
							lattice[partner_x+1][partner_y][partner_z] = my_object_array[t];
							lattice[partner_x][partner_y][partner_z] = null;
							// my_object_array_update
							my_object_array[r].xaxis = x+1;
							my_object_array[t].xaxis = partner_x+1;
						}
					}
				} else if(choice == my_object_array[r].bond_direction && x<(lattice_dim-2) && lattice[x+2][y][z] == null) {
					// bond direction is same as randomly chosen direction
					// All set to move
					// Lattice update
					lattice[partner_x+1][partner_y][partner_z] = lattice[partner_x][partner_y][partner_z];
					lattice[partner_x][partner_y][partner_z] = lattice[x][y][z];
					lattice[x][y][z] = null;
					// my_object_array update
					my_object_array[r].xaxis = x+1;
					my_object_array[t].xaxis = partner_x+1;
				}
			} else if (choice == 2) {
				// Left
				if (choice != my_object_array[r].bond_direction  && x>1 && lattice[x-1][y][z] == null) {
					if (my_object_array[r].bond_direction == 1) {
						// Same bond orientation in opposite direction
						// All set to move.
						// lattice update
						lattice[x-1][y][z] = my_object_array[r];
						lattice[x][y][z] = lattice[partner_x][partner_y][partner_z];
						lattice[partner_x][partner_y][partner_z] = null;
						// my_object_array update
						my_object_array[r].xaxis = x-1;
						my_object_array[t].xaxis = partner_x-1;
					} else {
						// Different bond orientation
						if (lattice[partner_x-1][partner_y][partner_z] == null) {
							// All set to move
							// lattice update
							lattice[x-1][y][z] = my_object_array[r];
							lattice[x][y][z] = null;
							lattice[partner_x-1][partner_y][partner_z] = my_object_array[t];
							lattice[partner_x][partner_y][partner_z] = null;
							// my_object_array_update
							my_object_array[r].xaxis = x-1;
							my_object_array[t].xaxis = partner_x-1;
						}
					}
				} else if(choice == my_object_array[r].bond_direction && x>2 && lattice[x-2][y][z] == null) {
					// bond direction is same as randomly chosen direction
					// All set to move
					// Lattice update
					lattice[partner_x-1][partner_y][partner_z] = lattice[partner_x][partner_y][partner_z];
					lattice[partner_x][partner_y][partner_z] = lattice[x][y][z];
					lattice[x][y][z] = null;
					// my_object_array update
					my_object_array[r].xaxis = x-1;
					my_object_array[t].xaxis = partner_x-1;
				}
			} else if (choice == 3) {
				// Up
				if (choice != my_object_array[r].bond_direction  && y<(lattice_dim-1) && lattice[x][y+1][z] == null) {
					if (my_object_array[r].bond_direction == 4) {
						// Same bond orientation in opposite direction
						// All set to move.
						// lattice update
						lattice[x][y+1][z] = my_object_array[r];
						lattice[x][y][z] = lattice[partner_x][partner_y][partner_z];
						lattice[partner_x][partner_y][partner_z] = null;
						// my_object_array update
						my_object_array[r].yaxis = y+1;
						my_object_array[t].yaxis = partner_y+1;
					} else {
						// Different bond orientation
						if (lattice[partner_x][partner_y+1][partner_z] == null) {
							// All set to move
							// lattice update
							lattice[x][y+1][z] = my_object_array[r];
							lattice[x][y][z] = null;
							lattice[partner_x][partner_y+1][partner_z] = my_object_array[t];
							lattice[partner_x][partner_y][partner_z] = null;
							// my_object_array_update
							my_object_array[r].yaxis = y+1;
							my_object_array[t].yaxis = partner_y+1;
						}
					}
				} else if(choice == my_object_array[r].bond_direction && y<(lattice_dim-2) && lattice[x][y+2][z] == null) {
					// bond direction is same as randomly chosen direction
					// All set to move
					// Lattice update
					lattice[partner_x][partner_y+1][partner_z] = lattice[partner_x][partner_y][partner_z];
					lattice[partner_x][partner_y][partner_z] = lattice[x][y][z];
					lattice[x][y][z] = null;
					// my_object_array update
					my_object_array[r].yaxis = y+1;
					my_object_array[t].yaxis = partner_y+1;
				}
			} else if (choice == 4) {
				// Down
				if (choice != my_object_array[r].bond_direction  && y>1 && lattice[x][y-1][z] == null) {
					if (my_object_array[r].bond_direction == 3) {
						// Same bond orientation in opposite direction
						// All set to move.
						// lattice update
						lattice[x][y-1][z] = my_object_array[r];
						lattice[x][y][z] = lattice[partner_x][partner_y][partner_z];
						lattice[partner_x][partner_y][partner_z] = null;
						// my_object_array update
						my_object_array[r].yaxis = y-1;
						my_object_array[t].yaxis = partner_y-1;
					} else {
						// Different bond orientation
						if (lattice[partner_x][partner_y-1][partner_z] == null) {
							// All set to move
							// lattice update
							lattice[x][y-1][z] = my_object_array[r];
							lattice[x][y][z] = null;
							lattice[partner_x][partner_y-1][partner_z] = my_object_array[t];
							lattice[partner_x][partner_y][partner_z] = null;
							// my_object_array_update
							my_object_array[r].yaxis = y-1;
							my_object_array[t].yaxis = partner_y-1;
						}
					}
				} else if(choice == my_object_array[r].bond_direction && y>2 && lattice[x][y-2][z] == null) {
					// bond direction is same as randomly chosen direction
					// All set to move
					// Lattice update
					lattice[partner_x][partner_y-1][partner_z] = lattice[partner_x][partner_y][partner_z];
					lattice[partner_x][partner_y][partner_z] = lattice[x][y][z];
					lattice[x][y][z] = null;
					// my_object_array update
					my_object_array[r].yaxis = y-1;
					my_object_array[t].yaxis = partner_y-1;
				}
			} else if (choice == 5) {
				// Front
				if (choice != my_object_array[r].bond_direction  && z<(lattice_dim-1) && lattice[x][y][z+1] == null) {
					if (my_object_array[r].bond_direction == 6) {
						// Same bond orientation in opposite direction
						// All set to move.
						// lattice update
						lattice[x][y][z+1] = my_object_array[r];
						lattice[x][y][z] = lattice[partner_x][partner_y][partner_z];
						lattice[partner_x][partner_y][partner_z] = null;
						// my_object_array update
						my_object_array[r].zaxis = z+1;
						my_object_array[t].zaxis = partner_z+1;
					} else {
						// Different bond orientation
						if (lattice[partner_x][partner_y][partner_z+1] == null) {
							// All set to move
							// lattice update
							lattice[x][y][z+1] = my_object_array[r];
							lattice[x][y][z] = null;
							lattice[partner_x][partner_y][partner_z+1] = my_object_array[t];
							lattice[partner_x][partner_y][partner_z] = null;
							// my_object_array_update
							my_object_array[r].zaxis = z+1;
							my_object_array[t].zaxis = partner_z+1;
						}
					}
				} else if(choice == my_object_array[r].bond_direction && z<(lattice_dim-2) && lattice[x][y][z+2] == null) {
					// bond direction is same as randomly chosen direction
					// All set to move
					// Lattice update
					lattice[partner_x][partner_y][partner_z+1] = lattice[partner_x][partner_y][partner_z];
					lattice[partner_x][partner_y][partner_z] = lattice[x][y][z];
					lattice[x][y][z] = null;
					// my_object_array update
					my_object_array[r].zaxis = z+1;
					my_object_array[t].zaxis = partner_z+1;
				}
			} else if (choice == 6) {
				// Back
				if (choice != my_object_array[r].bond_direction  && z>1 && lattice[x][y][z-1] == null) {
					if (my_object_array[r].bond_direction == 5) {
						// Same bond orientation in opposite direction
						// All set to move.
						// lattice update
						lattice[x][y][z-1] = my_object_array[r];
						lattice[x][y][z] = lattice[partner_x][partner_y][partner_z];
						lattice[partner_x][partner_y][partner_z] = null;
						// my_object_array update
						my_object_array[r].zaxis = z-1;
						my_object_array[t].zaxis = partner_z-1;
					} else {
						// Different bond orientation
						if (lattice[partner_x][partner_y][partner_z-1] == null) {
							// All set to move
							// lattice update
							lattice[x][y][z-1] = my_object_array[r];
							lattice[x][y][z] = null;
							lattice[partner_x][partner_y][partner_z-1] = my_object_array[t];
							lattice[partner_x][partner_y][partner_z] = null;
							// my_object_array_update
							my_object_array[r].zaxis = z-1;
							my_object_array[t].zaxis = partner_z-1;
						}
					}
				} else if(choice == my_object_array[r].bond_direction && z>2 && lattice[x][y][z-2] == null) {
					// bond direction is same as randomly chosen direction
					// All set to move
					// Lattice update
					lattice[partner_x][partner_y][partner_z-1] = lattice[partner_x][partner_y][partner_z];
					lattice[partner_x][partner_y][partner_z] = lattice[x][y][z];
					lattice[x][y][z] = null;
					// my_object_array update
					my_object_array[r].zaxis = z-1;
					my_object_array[t].zaxis = partner_z-1;
				}
			}
		}
	}
	
	public static void bond_formation_by_cytoc_apaf(int r)
	{
		int x = my_object_array[r].xaxis;
		int y = my_object_array[r].yaxis;
		int z = my_object_array[r].zaxis;
		if(my_object_array[r].molecule_bonded_with == 0)
		{
			if(Math.random() < p_apoptosome_form)
			{
				if(x != lattice_dim - 1)
				{
					if(lattice[x+1][y][z] != null)
					{
						//System.out.println(lattice[x][y][z]);
						//System.out.println(lattice[x+1][y][z]);
						if((lattice[x][y][z].molecule_type == 1  && lattice[x+1][y][z].molecule_type == 2) || (lattice[x][y][z].molecule_type == 2  && lattice[x+1][y][z].molecule_type == 1))
						{
							my_object_array[r].molecule_bonded_with = lattice[x+1][y][z].molecule_type;
	                        for (int j = 0; j < total_molecules; j++)
	                        {
	                        	if(my_object_array[j].xaxis == (x+1) && my_object_array[j].yaxis == y && my_object_array[z].zaxis == z)
	                        	{
	                        		my_object_array[j].molecule_bonded_with = lattice[x][y][z].molecule_type;
	                                //my_object_array[j].bond_direction = 2;
	                                //my_object_array[j].molecule_type = 3;
	                                break;
	                        	}
	                        }                            
	                    	lattice[x][y][z].molecule_bonded_with = lattice[x+1][y][z].molecule_type;
	                        lattice[x][y][z].bond_direction = 1;
	                        lattice[x+1][y][z].molecule_bonded_with = lattice[x][y][z].molecule_type;
	                        lattice[x+1][y][z].bond_direction = 2;
	                        lattice[x][y][z].molecule_type = 3;
	                        lattice[x+1][y][z].molecule_type = 3;
	                        my_object_array[r].bond_direction = 1;
	                        my_object_array[r].molecule_type = 3;
	                        apoptosome_count ++;                           
	                        for (int j = 0; j < total_molecules; j++)
	                        {
	                        	if(my_object_array[j].xaxis == (x+1) && my_object_array[j].yaxis == y && my_object_array[j].zaxis == z)
	                        	{
	                        		//my_object_array[j].molecule_bonded_with = lattice[x][y][z].molecule_type;
	                                my_object_array[j].bond_direction = 2;
	                                my_object_array[j].molecule_type = 3;
	                                break;
	                        	}
	                        }
	                        return;
						}
					}
				}
				if(y != lattice_dim - 1)
				{
					if(lattice[x][y+1][z] != null)
					{
						//System.out.println(lattice[x][y][z].molecule_type);
						//System.out.println(lattice[x][y+1][z].molecule_type);
						if((lattice[x][y][z].molecule_type == 1  && lattice[x][y+1][z].molecule_type == 2) || (lattice[x][y][z].molecule_type == 2  && lattice[x][y+1][z].molecule_type == 1))
						{
							my_object_array[r].molecule_bonded_with = lattice[x][y+1][z].molecule_type;
	                        for (int j = 0; j < total_molecules; j++)
	                        {
	                        	if(my_object_array[j].xaxis == x && my_object_array[j].yaxis == (y+1) && my_object_array[j].zaxis == z)
	                        	{
	                        		my_object_array[j].molecule_bonded_with = lattice[x][y][z].molecule_type;
	                                //my_object_array[j].bond_direction = 4;
	                                //my_object_array[j].molecule_type = 3;
	                                break;
	                        	}
	                        }
	                    	lattice[x][y][z].molecule_bonded_with = lattice[x][y+1][z].molecule_type;
	                        lattice[x][y][z].bond_direction = 3;
	                        lattice[x][y+1][z].molecule_bonded_with = lattice[x][y][z].molecule_type;
	                        lattice[x][y+1][z].bond_direction = 4;
	                        lattice[x][y][z].molecule_type = 3;
	                        lattice[x][y+1][z].molecule_type = 3;
	                        my_object_array[r].bond_direction = 3;
	                        my_object_array[r].molecule_type = 3;
	                        apoptosome_count ++;
	                        for (int j = 0; j < total_molecules; j++)
	                        {
	                        	if(my_object_array[j].xaxis == x && my_object_array[j].yaxis == (y+1) && my_object_array[j].zaxis == z)
	                        	{
	                        		//my_object_array[j].molecule_bonded_with = lattice[x][y][z].molecule_type;
	                                my_object_array[j].bond_direction = 4;
	                                my_object_array[j].molecule_type = 3;
	                                break;
	                        	}
	                        }
	                        return;
						}
					}
				}
				if(z != lattice_dim - 1)
				{
					if(lattice[x][y][z+1] != null)
					{
						if((lattice[x][y][z].molecule_type == 1  && lattice[x][y][z+1].molecule_type == 2) || (lattice[x][y][z].molecule_type == 2  && lattice[x][y][z+1].molecule_type == 1))
						{
							my_object_array[r].molecule_bonded_with = lattice[x][y][z+1].molecule_type;
	                        for (int j = 0; j < total_molecules; j++)
	                        {
	                        	if(my_object_array[j].xaxis == x && my_object_array[j].yaxis == y && my_object_array[j].zaxis == z+1)
	                        	{
	                        		my_object_array[j].molecule_bonded_with = lattice[x][y][z].molecule_type;
	                                //my_object_array[j].bond_direction = 6;
	                                //my_object_array[j].molecule_type = 3;
	                                break;
	                        	}
	                        }
	                    	lattice[x][y][z].molecule_bonded_with = lattice[x][y][z+1].molecule_type;
	                        lattice[x][y][z].bond_direction = 5;
	                        lattice[x][y][z+1].molecule_bonded_with = lattice[x][y][z].molecule_type;
	                        lattice[x][y][z+1].bond_direction = 6;
	                        lattice[x][y][z].molecule_type = 3;
	                        lattice[x][y][z+1].molecule_type = 3;
	                        my_object_array[r].bond_direction = 5;
	                        my_object_array[r].molecule_type = 3;
	                        apoptosome_count ++;
	                        for (int j = 0; j < total_molecules; j++)
	                        {
	                        	if(my_object_array[j].xaxis == x && my_object_array[j].yaxis == y && my_object_array[j].zaxis == z+1)
	                        	{
	                        		//my_object_array[j].molecule_bonded_with = lattice[x][y][z].molecule_type;
	                                my_object_array[j].bond_direction = 6;
	                                my_object_array[j].molecule_type = 3;
	                                break;
	                        	}
	                        }
	                        return;
						}
					}
				}
	            if (x > 0)
	            {
	            	if(lattice[x-1][y][z] != null)
	            	{
	            		if ((lattice[x][y][z].molecule_type == 1  && lattice[x-1][y][z].molecule_type == 2) || (lattice[x][y][z].molecule_type == 2  && lattice[x-1][y][z].molecule_type == 1))
	            		{
	                        my_object_array[r].molecule_bonded_with = lattice[x-1][y][z].molecule_type;
	                        for (int j = 0; j < total_molecules; j++)
	                        {
	                        	if(my_object_array[j].xaxis == x-1 & my_object_array[j].yaxis == y & my_object_array[j].zaxis == z)
	                        	{
	                                my_object_array[j].molecule_bonded_with = lattice[x][y][z].molecule_type;
	                                //my_object_array[j].bond_direction = 1;
	                                //my_object_array[j].molecule_type = 3;
	                                break;
	                        	}
	                        }                     
	                        lattice[x][y][z].molecule_bonded_with = lattice[x - 1][y][z].molecule_type;
	                        lattice[x][y][z].bond_direction = 2;
	                        lattice[x - 1][y][z].molecule_bonded_with = lattice[x][y][z].molecule_type;
	                        lattice[x - 1][y][z].bond_direction = 1;
	                        lattice[x][y][z].molecule_type = 3;
	                        lattice[x - 1][y][z].molecule_type = 3;
	                        my_object_array[r].bond_direction = 2;
	                        my_object_array[r].molecule_type = 3;
	                        apoptosome_count ++;
	                        for (int j = 0; j < total_molecules; j++)
	                        {
	                        	if(my_object_array[j].xaxis == x-1 & my_object_array[j].yaxis == y & my_object_array[j].zaxis == z)
	                        	{
	                                //my_object_array[j].molecule_bonded_with = lattice[x][y][z].molecule_type;
	                                my_object_array[j].bond_direction = 1;
	                                my_object_array[j].molecule_type = 3;
	                                break;
	                        	}
	                        } 
	                        return;
	            		}
	            	}	
	            }
	            if (y > 0)
	            {
	            	if(lattice[x][y-1][z] != null)
	            	{
	            		if ((lattice[x][y][z].molecule_type == 1  && lattice[x][y-1][z].molecule_type == 2) || (lattice[x][y][z].molecule_type == 2  && lattice[x][y-1][z].molecule_type == 1))
	            		{
	                        my_object_array[r].molecule_bonded_with = lattice[x][y-1][z].molecule_type;
	                        for (int j = 0; j < total_molecules; j++)
	                        {
	                        	if(my_object_array[j].xaxis == x & my_object_array[j].yaxis == y-1 & my_object_array[j].zaxis == z)
	                        	{
	                                my_object_array[j].molecule_bonded_with = lattice[x][y][z].molecule_type;
	                                //my_object_array[j].bond_direction = 3;
	                                //my_object_array[j].molecule_type = 3;
	                                break;
	                        	}
	                        }                     
	                        lattice[x][y][z].molecule_bonded_with = lattice[x][y-1][z].molecule_type;
	                        lattice[x][y][z].bond_direction = 4;
	                        lattice[x][y-1][z].molecule_bonded_with = lattice[x][y][z].molecule_type;
	                        lattice[x][y-1][z].bond_direction = 3;
	                        lattice[x][y][z].molecule_type = 3;
	                        lattice[x][y-1][z].molecule_type = 3;
	                        my_object_array[r].bond_direction = 4;
	                        my_object_array[r].molecule_type = 3;
	                        apoptosome_count ++;
	                        for (int j = 0; j < total_molecules; j++)
	                        {
	                        	if(my_object_array[j].xaxis == x & my_object_array[j].yaxis == y-1 & my_object_array[j].zaxis == z)
	                        	{
	                                //my_object_array[j].molecule_bonded_with = lattice[x][y][z].molecule_type;
	                                my_object_array[j].bond_direction = 3;
	                                my_object_array[j].molecule_type = 3;
	                                break;
	                        	}
	                        } 
	                        return;
	            		}
	            	}	
	            }	
	            if (z > 0)
	            {
	            	if(lattice[x][y][z-1] != null)
	            	{
	            		if ((lattice[x][y][z].molecule_type == 1  && lattice[x][y][z-1].molecule_type == 2) || (lattice[x][y][z].molecule_type == 2  && lattice[x][y][z-1].molecule_type == 1))
	            		{
	                        my_object_array[r].molecule_bonded_with = lattice[x][y][z-1].molecule_type;
	                        for (int j = 0; j < total_molecules; j++)
	                        {
	                        	if(my_object_array[j].xaxis == x & my_object_array[j].yaxis == y & my_object_array[j].zaxis == z-1)
	                        	{
	                                my_object_array[j].molecule_bonded_with = lattice[x][y][z].molecule_type;
	                                //my_object_array[j].bond_direction = 5;
	                                //my_object_array[j].molecule_type = 3;
	                                break;
	                        	}
	                        }                     
	                        lattice[x][y][z].molecule_bonded_with = lattice[x][y][z-1].molecule_type;
	                        lattice[x][y][z].bond_direction = 6;
	                        lattice[x][y][z-1].molecule_bonded_with = lattice[x][y][z].molecule_type;
	                        lattice[x][y][z-1].bond_direction = 5;
	                        lattice[x][y][z].molecule_type = 3;
	                        lattice[x][y][z-1].molecule_type = 3;
	                        my_object_array[r].bond_direction = 6;
	                        my_object_array[r].molecule_type = 3;
	                        apoptosome_count ++;
	                        for (int j = 0; j < total_molecules; j++)
	                        {
	                        	if(my_object_array[j].xaxis == x & my_object_array[j].yaxis == y & my_object_array[j].zaxis == z-1)
	                        	{
	                                //my_object_array[j].molecule_bonded_with = lattice[x][y][z].molecule_type;
	                                my_object_array[j].bond_direction = 5;
	                                my_object_array[j].molecule_type = 3;
	                                break;
	                        	}
	                        }  
	                        return;
	            		}
	            	}	
	            }
			}
		}
	}

	public static void bond_formation_by_apoptosome(int r)
	{
		
  	  int x=my_object_array[r].xaxis;
  	  int y=my_object_array[r].yaxis;
  	  int z=my_object_array[r].zaxis;
  	  
  	  if(my_object_array[r].molecule_type==3)
  	  {
  	 int  direction_to_be_excluded = my_object_array[r].bond_direction;
  	 
  	 for(int j=1;j<6;j++)
  	 {
  	 if(j!=direction_to_be_excluded)
  	 
  	 
  	 if(j==1)
  	 {	 
  		 
  	 if(0<x+1&&x+1 <= (lattice_dim - 1)&&lattice[x+1][y][z]!=null && lattice[x+1][y][z].molecule_type==3)
  	 {
  		 
  	 lattice[x][y][z].molecule_type=4;
  	 my_object_array[r].molecule_type=4;
  	 lattice[x+1][y][z].molecule_type=4;
  	 dimer_count++;
  	 for (int k = 0;k<total_molecules;k++)
  	 {
                               if(my_object_array[k].xaxis == x+1 && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z)
                                   my_object_array[k].molecule_type=4;
  	 }
  	 
  	 if(direction_to_be_excluded==2)
  	 {	 
  		
  	     lattice[x-1][y][z].molecule_type=4;
  	 
  	 for (int k = 0;k<total_molecules;k++)
	 {
                           if(my_object_array[k].xaxis == x-1 && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z)
                               my_object_array[k].molecule_type=4;
	 }
	 
  	 }
  	 else if(direction_to_be_excluded==3)
  	 {	 
  		
      	     lattice[x][y+1][z].molecule_type=4;
      	 
      	 for (int k = 0;k<total_molecules;k++)
  	 {
                               if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y+1 && my_object_array[k].zaxis == z)
                                   my_object_array[k].molecule_type=4;
  	 }
  	 
      	 }
  	 else if(direction_to_be_excluded==4)
  	 {	 
  		
      	     lattice[x][y-1][z].molecule_type=4;
      	 
      	 for (int k = 0;k<total_molecules;k++)
  	 {
                               if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y-1 && my_object_array[k].zaxis == z)
                                   my_object_array[k].molecule_type=4;
  	 }
  	 
      	 }
  	 else if(direction_to_be_excluded==5)
  	 {	 
  		
      	     lattice[x][y][z+1].molecule_type=4;
      	 
      	 for (int k = 0;k<total_molecules;k++)
  	 {
                               if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z+1)
                                   my_object_array[k].molecule_type=4;
  	 }
  	 
      	 }
  	 else if(direction_to_be_excluded==6)
  	 {	 
  		
      	     lattice[x][y][z-1].molecule_type=4;
      	 
      	 for (int k = 0;k<total_molecules;k++)
  	 {
                               if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z-1)
                                   my_object_array[k].molecule_type=4;
  	 }
  	 
      	 }
  	
  	int direction_of_other_apoptosome=lattice[x+1][y][z].bond_direction;
  	
  	if(direction_of_other_apoptosome==1)
 	 {	 
 	     lattice[x+2][y][z].molecule_type=4;
 	 
 	 for (int k = 0;k<total_molecules;k++)
	 {
                          if(my_object_array[k].xaxis == x+2 && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z)
                              my_object_array[k].molecule_type=4;
	 }
	 
 	 }
  	
 	 else if(direction_of_other_apoptosome==3)
 	 {	 
     	     lattice[x+1][y+1][z].molecule_type=4;
     	 
     	 for (int k = 0;k<total_molecules;k++)
 	 {
                              if(my_object_array[k].xaxis == x+1 && my_object_array[k].yaxis == y+1 && my_object_array[k].zaxis == z)
                                  my_object_array[k].molecule_type=4;
 	 }
 	 
     	 }
 	 else if(direction_of_other_apoptosome==4)
 	 {	 
     	     lattice[x+1][y-1][z].molecule_type=4;
     	 
     	 for (int k = 0;k<total_molecules;k++)
 	 {
                              if(my_object_array[k].xaxis == x+1 && my_object_array[k].yaxis == y-1 && my_object_array[k].zaxis == z)
                                  my_object_array[k].molecule_type=4;
 	 }
 	 
     	 }
 	 else if(direction_of_other_apoptosome==5)
 	 {	 
     	     lattice[x+1][y][z+1].molecule_type=4;
     	 
     	 for (int k = 0;k<total_molecules;k++)
 	 {
                              if(my_object_array[k].xaxis == x+1 && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z+1)
                                  my_object_array[k].molecule_type=4;
 	 }
 	 
     	 }
 	 else if(direction_of_other_apoptosome==6)
 	 {	 
     	     lattice[x+1][y][z-1].molecule_type=4;
     	 
     	 for (int k = 0;k<total_molecules;k++)
 	 {
                              if(my_object_array[k].xaxis == x+1 && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z-1)
                                  my_object_array[k].molecule_type=4;
 	 }
 	 
     	 }
  	
  	
  	 
  	 }
  	 break;
  	  }
  	 
  	 
  	 
  	 else if (j==2)
  	        {
  	
  	    
	 if(0<=x-1&&x-1 <= (lattice_dim - 1) && lattice[x-1][y][z]!=null && lattice[x-1][y][z].molecule_type==3)
	 {
	 System.out.println("Bumba");		 
	 lattice[x][y][z].molecule_type=4;
	 my_object_array[r].molecule_type=4;
	 lattice[x-1][y][z].molecule_type=4;
	 dimer_count++;
	 for (int k = 0;k<total_molecules;k++)
	 {
                          if(my_object_array[k].xaxis == x-1 && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z)
                              my_object_array[k].molecule_type=4;
	 }
	 
	 if(direction_to_be_excluded==1)
	 {	 
		 
	     lattice[x+1][y][z].molecule_type=4;
	 
	 for (int k = 0;k<total_molecules;k++)
	 {
                      if(my_object_array[k].xaxis == x+1 && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z)
                          my_object_array[k].molecule_type=4;
	 }
	 
	 }
	 else if(direction_to_be_excluded==3)
	 {	 
 	     lattice[x][y+1][z].molecule_type=4;
 	 
 	 for (int k = 0;k<total_molecules;k++)
	 {
                          if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y+1 && my_object_array[k].zaxis == z)
                              my_object_array[k].molecule_type=4;
	 }
	 
 	 }
	 else if(direction_to_be_excluded==4)
	 {	 
 	     lattice[x][y-1][z].molecule_type=4;
 	 
 	 for (int k = 0;k<total_molecules;k++)
	 {
                          if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y-1 && my_object_array[k].zaxis == z)
                              my_object_array[k].molecule_type=4;
	 }
	 
 	 }
	 else if(direction_to_be_excluded==5)
	 {	 
 	     lattice[x][y][z+1].molecule_type=4;
 	 
 	 for (int k = 0;k<total_molecules;k++)
	 {
                          if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z+1)
                              my_object_array[k].molecule_type=4;
	 }
	 
 	 }
	 else if(direction_to_be_excluded==6)
	 {	 
 	     lattice[x][y][z-1].molecule_type=4;
 	 
 	 for (int k = 0;k<total_molecules;k++)
	 {
                          if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z-1)
                              my_object_array[k].molecule_type=4;
	 }
	 
 	 }
	 
	int direction_of_other_apoptosome=lattice[x-1][y][z].bond_direction;
	
	if(direction_of_other_apoptosome==2)
	 {	 
	     lattice[x-2][y][z].molecule_type=4;
	 
	 for (int k = 0;k<total_molecules;k++)
	 {
                     if(my_object_array[k].xaxis == x-2 && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z)
                         my_object_array[k].molecule_type=4;
	 }
	 
	 }
	
	 else if(direction_of_other_apoptosome==3)
	 {	 
	     lattice[x-1][y+1][z].molecule_type=4;
	 
	 for (int k = 0;k<total_molecules;k++)
	 {
                         if(my_object_array[k].xaxis == x-1 && my_object_array[k].yaxis == y+1 && my_object_array[k].zaxis == z)
                             my_object_array[k].molecule_type=4;
	 }
	 
	 }
	 else if(direction_of_other_apoptosome==4)
	 {	 
	     lattice[x-1][y-1][z].molecule_type=4;
	 
	 for (int k = 0;k<total_molecules;k++)
	 {
                         if(my_object_array[k].xaxis == x-1 && my_object_array[k].yaxis == y-1 && my_object_array[k].zaxis == z)
                             my_object_array[k].molecule_type=4;
	 }
	 
	 }
	 else if(direction_of_other_apoptosome==5)
	 {	 
	     lattice[x-1][y][z+1].molecule_type=4;
	 
	 for (int k = 0;k<total_molecules;k++)
	 {
                         if(my_object_array[k].xaxis == x-1 && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z+1)
                             my_object_array[k].molecule_type=4;
	 }
	 
	 }
	 else if(direction_of_other_apoptosome==6)
	 {	 
	     lattice[x-1][y][z-1].molecule_type=4;
	 
	 for (int k = 0;k<total_molecules;k++)
	 {
                         if(my_object_array[k].xaxis == x-1 && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z-1)
                             my_object_array[k].molecule_type=4;
	 }
	 
	 }
	 }
	 break;
  	        }
	 
	 
	 
	 else if (j==3)
	    	        {
	    	
	    	    	 
	 if(0<=y+1&&y+1 != (lattice_dim - 1)&&lattice[x][y+1][z]!=null && lattice[x][y+1][z].molecule_type==3)
	 {
	 lattice[x][y][z].molecule_type=4;
	 my_object_array[r].molecule_type=4;
	 lattice[x][y+1][z].molecule_type=4;
	 dimer_count++;
	 for (int k = 0;k<total_molecules;k++)
	 {
	                            if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y+1 && my_object_array[k].zaxis == z)
	                                my_object_array[k].molecule_type=4;
	 }
	 
	 if(direction_to_be_excluded==1)
	 {	 
	     lattice[x+1][y][z].molecule_type=4;
	 
	 for (int k = 0;k<total_molecules;k++)
	 {
	                        if(my_object_array[k].xaxis == x+1 && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z)
	                            my_object_array[k].molecule_type=4;
	 }
	 
	 }
	 else if(direction_to_be_excluded==2)
	 {	 
	   	     lattice[x-1][y][z].molecule_type=4;
	   	 
	   	 for (int k = 0;k<total_molecules;k++)
	 {
	                            if(my_object_array[k].xaxis == x-1 && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z)
	                                my_object_array[k].molecule_type=4;
	 }
	 
	   	 }
	 else if(direction_to_be_excluded==4)
	 {	 
	   	     lattice[x][y-1][z].molecule_type=4;
	   	 
	   	 for (int k = 0;k<total_molecules;k++)
	 {
	                            if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y-1 && my_object_array[k].zaxis == z)
	                                my_object_array[k].molecule_type=4;
	 }
	 
	   	 }
	 else if(direction_to_be_excluded==5)
	 {	 
	   	     lattice[x][y][z+1].molecule_type=4;
	   	 
	   	 for (int k = 0;k<total_molecules;k++)
	 {
	                            if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z+1)
	                                my_object_array[k].molecule_type=4;
	 }
	 
	   	 }
	 else if(direction_to_be_excluded==6)
	 {	 
	   	     lattice[x][y][z-1].molecule_type=4;
	   	 
	   	 for (int k = 0;k<total_molecules;k++)
	 {
	                            if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z-1)
	                                my_object_array[k].molecule_type=4;
	 }
	 
	   	 }
	 
	int direction_of_other_apoptosome=lattice[x][y+1][z].bond_direction;
	
	
	if(direction_of_other_apoptosome==1)
	 {	 
	     lattice[x+1][y+1][z].molecule_type=4;
	 
	 for (int k = 0;k<total_molecules;k++)
	 {
	                       if(my_object_array[k].xaxis == x+1 && my_object_array[k].yaxis == y+1 && my_object_array[k].zaxis == z)
	                           my_object_array[k].molecule_type=4;
	 }
	 
	 }
	
	if(direction_of_other_apoptosome==2)
	 {	 
	     lattice[x-1][y+1][z].molecule_type=4;
	 
	 for (int k = 0;k<total_molecules;k++)
	 {
	                       if(my_object_array[k].xaxis == x-1 && my_object_array[k].yaxis == y+1 && my_object_array[k].zaxis == z)
	                           my_object_array[k].molecule_type=4;
	 }
	 
	 }
	
	 else if(direction_of_other_apoptosome==3)
	 {	 
	  	     lattice[x][y+2][z].molecule_type=4;
	  	 
	  	 for (int k = 0;k<total_molecules;k++)
	 {
	                           if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y+2 && my_object_array[k].zaxis == z)
	                               my_object_array[k].molecule_type=4;
	 }
	 
	  	 }
	 
	 else if(direction_of_other_apoptosome==5)
	 {	 
	  	     lattice[x][y+1][z+1].molecule_type=4;
	  	 
	  	 for (int k = 0;k<total_molecules;k++)
	 {
	                           if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y+1 && my_object_array[k].zaxis == z+1)
	                               my_object_array[k].molecule_type=4;
	 }
	 
	  	 }
	 else if(direction_of_other_apoptosome==6)
	 {	 
	  	     lattice[x][y+1][z-1].molecule_type=4;
	  	 
	  	 for (int k = 0;k<total_molecules;k++)
	 {
	                           if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y+1 && my_object_array[k].zaxis == z-1)
	                               my_object_array[k].molecule_type=4;
	 }
	 
	  	 }
	 }	 
	 
	 break;
	    	        }
	 else if (j==4)
	    	        {
	    	
	    	    	 
	 if(0<=y-1&&y-1 != (lattice_dim - 1)&& lattice[x][y-1][z]!=null && lattice[x][y-1][z].molecule_type==3)
	 {
	 lattice[x][y][z].molecule_type=4;
	 my_object_array[r].molecule_type=4;
	 lattice[x][y-1][z].molecule_type=4;
	 dimer_count++;
	 for (int k = 0;k<total_molecules;k++)
	 {
	                            if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y-1 && my_object_array[k].zaxis == z)
	                                my_object_array[k].molecule_type=4;
	 }
	 
	 if(direction_to_be_excluded==1)
	 {	 
	     lattice[x+1][y][z].molecule_type=4;
	 
	 for (int k = 0;k<total_molecules;k++)
	 {
	                        if(my_object_array[k].xaxis == x+1 && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z)
	                            my_object_array[k].molecule_type=4;
	 }
	 
	 }
	 else if(direction_to_be_excluded==2)
	 {	 
	   	     lattice[x-1][y][z].molecule_type=4;
	   	 
	   	 for (int k = 0;k<total_molecules;k++)
	 {
	                            if(my_object_array[k].xaxis == x-1 && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z)
	                                my_object_array[k].molecule_type=4;
	 }
	 
	   	 }
	 else if(direction_to_be_excluded==3)
	 {	 
	   	     lattice[x][y+1][z].molecule_type=4;
	   	 
	   	 for (int k = 0;k<total_molecules;k++)
	 {
	                            if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y+1 && my_object_array[k].zaxis == z)
	                                my_object_array[k].molecule_type=4;
	 }
	 
	   	 }
	 else if(direction_to_be_excluded==5)
	 {	 
	   	     lattice[x][y][z+1].molecule_type=4;
	   	 
	   	 for (int k = 0;k<total_molecules;k++)
	 {
	                            if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z+1)
	                                my_object_array[k].molecule_type=4;
	 }
	 
	   	 }
	 else if(direction_to_be_excluded==6)
	 {	 
	   	     lattice[x][y][z-1].molecule_type=4;
	   	 
	   	 for (int k = 0;k<total_molecules;k++)
	 {
	                            if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z-1)
	                                my_object_array[k].molecule_type=4;
	 }
	 
	   	 }
	 
	int direction_of_other_apoptosome=lattice[x][y-1][z].bond_direction;
	
	
	if(direction_of_other_apoptosome==1)
	 {	 
	     lattice[x+1][y-1][z].molecule_type=4;
	 
	 for (int k = 0;k<total_molecules;k++)
	 {
	                       if(my_object_array[k].xaxis == x+1 && my_object_array[k].yaxis == y-1 && my_object_array[k].zaxis == z)
	                           my_object_array[k].molecule_type=4;
	 }
	 
	 }
	
	if(direction_of_other_apoptosome==2)
	 {	 
	     lattice[x-1][y-1][z].molecule_type=4;
	 
	 for (int k = 0;k<total_molecules;k++)
	 {
	                       if(my_object_array[k].xaxis == x-1 && my_object_array[k].yaxis == y-1 && my_object_array[k].zaxis == z)
	                           my_object_array[k].molecule_type=4;
	 }
	 
	 }
	
	 else if(direction_of_other_apoptosome==4)
	 {	 
	  	     lattice[x][y-2][z].molecule_type=4;
	  	 
	  	 for (int k = 0;k<total_molecules;k++)
	 {
	                           if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y-2 && my_object_array[k].zaxis == z)
	                               my_object_array[k].molecule_type=4;
	 }
	 
	  	 }
	 
	 else if(direction_of_other_apoptosome==5)
	 {	 
	  	     lattice[x][y-1][z+1].molecule_type=4;
	  	 
	  	 for (int k = 0;k<total_molecules;k++)
	 {
	                           if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y-1 && my_object_array[k].zaxis == z+1)
	                               my_object_array[k].molecule_type=4;
	 }
	 
	  	 }
	 else if(direction_of_other_apoptosome==6)
	 {	 
	  	     lattice[x][y-1][z-1].molecule_type=4;
	  	 
	  	 for (int k = 0;k<total_molecules;k++)
	 {
	                           if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y-1 && my_object_array[k].zaxis == z-1)
	                               my_object_array[k].molecule_type=4;
	 }
	 
	  	 }
	 }	 
	 
	 break;
	    	        } 
	 
	 else if(j==5)
	 {
	    	
	    	    	 
	 if(0<=z+1&&z+1 != (lattice_dim - 1)&&lattice[x][y][z+1]!=null && lattice[x][y][z+1].molecule_type==3)
	 {
	 lattice[x][y][z].molecule_type=4;
	 my_object_array[r].molecule_type=4;
	 lattice[x][y][z+1].molecule_type=4;
	 dimer_count++;
	 for (int k = 0;k<total_molecules;k++)
	 {
	                            if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z+1)
	                                my_object_array[k].molecule_type=4;
	 }
	 
	 if(direction_to_be_excluded==1)
	 {	 
	     lattice[x+1][y][z].molecule_type=4;
	 
	 for (int k = 0;k<total_molecules;k++)
	 {
	                        if(my_object_array[k].xaxis == x+1 && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z)
	                            my_object_array[k].molecule_type=4;
	 }
	 
	 }
	 else if(direction_to_be_excluded==2)
	 {	 
	   	     lattice[x-1][y][z].molecule_type=4;
	   	 
	   	 for (int k = 0;k<total_molecules;k++)
	 {
	                            if(my_object_array[k].xaxis == x-1 && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z)
	                                my_object_array[k].molecule_type=4;
	 }
	 
	   	 }
	 else if(direction_to_be_excluded==3)
	 {	 
	   	     lattice[x][y+1][z].molecule_type=4;
	   	 
	   	 for (int k = 0;k<total_molecules;k++)
	 {
	                            if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y+1 && my_object_array[k].zaxis == z)
	                                my_object_array[k].molecule_type=4;
	 }
	 
	   	 }
	 else if(direction_to_be_excluded==4)
	 {	 
	   	     lattice[x][y-1][z].molecule_type=4;
	   	 
	   	 for (int k = 0;k<total_molecules;k++)
	 {
	                            if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y-1 && my_object_array[k].zaxis == z)
	                                my_object_array[k].molecule_type=4;
	 }
	 
	   	 }
	 else if(direction_to_be_excluded==6)
	 {	 
	   	  
		 if(lattice[x][y][z-1]!=null && z!=lattice_dim-1)
				 lattice[x][y][z-1].molecule_type=4;
	   	 
	   	 for (int k = 0;k<total_molecules;k++)
	 {
	                            if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z-1)
	                                my_object_array[k].molecule_type=4;
	 }
	 
	   	 }
	 
	int direction_of_other_apoptosome=lattice[x][y][z+1].bond_direction;
	
	
	if(direction_of_other_apoptosome==1)
	 {	 
	     lattice[x+1][y][z+1].molecule_type=4;
	 
	 for (int k = 0;k<total_molecules;k++)
	 {
	                       if(my_object_array[k].xaxis == x+1 && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z+1)
	                           my_object_array[k].molecule_type=4;
	 }
	 
	 }
	
	else if(direction_of_other_apoptosome==2)
	 {	 
	     lattice[x-1][y][z+1].molecule_type=4;
	 
	 for (int k = 0;k<total_molecules;k++)
	 {
	                       if(my_object_array[k].xaxis == x-1 && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z+1)
	                           my_object_array[k].molecule_type=4;
	 }
	 
	 }
	else if(direction_of_other_apoptosome==3)
	 {	 
	     lattice[x][y+1][z+1].molecule_type=4;
	 
	 for (int k = 0;k<total_molecules;k++)
	 {
	                       if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y+1 && my_object_array[k].zaxis == z-1)
	                           my_object_array[k].molecule_type=4;
	 }
	 
	 }
	 else if(direction_of_other_apoptosome==4)
	 {	 
	  	     lattice[x][y-1][z+1].molecule_type=4;
	  	 
	  	 for (int k = 0;k<total_molecules;k++)
	 {
	                           if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y-1 && my_object_array[k].zaxis == z+1)
	                               my_object_array[k].molecule_type=4;
	 }
	 
	  	 }
	 
	 else if(direction_of_other_apoptosome==5)
	 {	 
	  	     lattice[x][y][z+2].molecule_type=4;
	  	 
	  	 for (int k = 0;k<total_molecules;k++)
	 {
	                           if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z+2)
	                               my_object_array[k].molecule_type=4;
	 }
	 
	  	 }
	 
	 }
	 break;
	 } 
	 
	 
	 else if(j==6)
	 {
	    	
	    	    	 
	 if(0<=z-1&&z-1 != (lattice_dim - 1)&&lattice[x][y][z-1]!=null && lattice[x][y][z-1].molecule_type==3)
	 {
	 lattice[x][y][z].molecule_type=4;
	 my_object_array[r].molecule_type=4;
	 lattice[x][y][z-1].molecule_type=4;
	 dimer_count++;
	 for (int k = 0;k<total_molecules;k++)
	 {
	                            if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z-1)
	                                my_object_array[k].molecule_type=4;
	 }
	 
	 if(direction_to_be_excluded==1)
	 {	 
	     lattice[x+1][y][z].molecule_type=4;
	 
	 for (int k = 0;k<total_molecules;k++)
	 {
	                        if(my_object_array[k].xaxis == x+1 && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z)
	                            my_object_array[k].molecule_type=4;
	 }
	 
	 }
	 else if(direction_to_be_excluded==2)
	 {	 
	   	     lattice[x-1][y][z].molecule_type=4;
	   	 
	   	 for (int k = 0;k<total_molecules;k++)
	 {
	                            if(my_object_array[k].xaxis == x-1 && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z)
	                                my_object_array[k].molecule_type=4;
	 }
	 
	   	 }
	 else if(direction_to_be_excluded==3)
	 {	 
	   	     lattice[x][y+1][z].molecule_type=4;
	   	 
	   	 for (int k = 0;k<total_molecules;k++)
	 {
	                            if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y+1 && my_object_array[k].zaxis == z)
	                                my_object_array[k].molecule_type=4;
	 }
	 
	   	 }
	 else if(direction_to_be_excluded==4)
	 {	 
	   	     lattice[x][y-1][z].molecule_type=4;
	   	 
	   	 for (int k = 0;k<total_molecules;k++)
	 {
	                            if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y-1 && my_object_array[k].zaxis == z)
	                                my_object_array[k].molecule_type=4;
	 }
	 
	   	 }
	 else if(direction_to_be_excluded==5)
	 {	 
	   	     lattice[x][y][z+1].molecule_type=4;
	   	 
	   	 for (int k = 0;k<total_molecules;k++)
	 {
	                            if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z+1)
	                                my_object_array[k].molecule_type=4;
	 }
	 
	   	 }
	 
	int direction_of_other_apoptosome=lattice[x][y][z-1].bond_direction;
	
	
	if(direction_of_other_apoptosome==1)
	 {	 
	     lattice[x+1][y][z-1].molecule_type=4;
	 
	 for (int k = 0;k<total_molecules;k++)
	 {
	                       if(my_object_array[k].xaxis == x+1 && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z-1)
	                           my_object_array[k].molecule_type=4;
	 }
	 
	 }
	
	else if(direction_of_other_apoptosome==2)
	 {	 
	     lattice[x-1][y][z-1].molecule_type=4;
	 
	 for (int k = 0;k<total_molecules;k++)
	 {
	                       if(my_object_array[k].xaxis == x-1 && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z-1)
	                           my_object_array[k].molecule_type=4;
	 }
	 
	 }
	else if(direction_of_other_apoptosome==3)
	 {	 
	     lattice[x][y+1][z-1].molecule_type=4;
	 
	 for (int k = 0;k<total_molecules;k++)
	 {
	                       if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y+1 && my_object_array[k].zaxis == z-1)
	                           my_object_array[k].molecule_type=4;
	 }
	 
	 }
	 else if(direction_of_other_apoptosome==4)
	 {	 
	  	     lattice[x][y-1][z-1].molecule_type=4;
	  	 
	  	 for (int k = 0;k<total_molecules;k++)
	 {
	                           if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y-1 && my_object_array[k].zaxis == z-1)
	                               my_object_array[k].molecule_type=4;
	 }
	 
	  	 }
	 
	 else if(direction_of_other_apoptosome==6)
	 {	 
	  	     lattice[x][y][z+2].molecule_type=4;
	  	 
	  	 for (int k = 0;k<total_molecules;k++)
	 {
	                           if(my_object_array[k].xaxis == x && my_object_array[k].yaxis == y && my_object_array[k].zaxis == z+2)
	                               my_object_array[k].molecule_type=4;
	 }
	 
	  	 }
	 
	
	
	 }
	 break;
	 }
  	 
	 }
	                       
	 
		 
	 }
  	 
  	 
  	  }
	
	public static void dissociation_apoptosome(int r)
	{
		int x = my_object_array[r].xaxis;
		int y = my_object_array[r].yaxis;
		int z = my_object_array[r].zaxis;
		if(my_object_array[r].molecule_type == 3)
		{
			
			//Random rand_double = new Random();
			if((Math.random() < p_apoptosome_break) && (apoptosome_count > 0))
			{
				// System.out.println(lattice[x+1][y][z]);
				// System.out.println(lattice[x][y][z]);
				if(lattice[x][y][z].bond_direction == 1 && x != (lattice_dim-1) && lattice[x+1][y][z] != null)
				{
					lattice[x][y][z].molecule_type = lattice[x+1][y][z].molecule_bonded_with;
					for(int j = 0; j < total_molecules; j++)
					{
						if(my_object_array[j].xaxis == x && my_object_array[j].yaxis == y && my_object_array[j].zaxis == z)
						{
							my_object_array[j].molecule_type = lattice[x+1][y][z].molecule_bonded_with;
                            break;	
						}
					}
					lattice[x+1][y][z].molecule_type = lattice[x][y][z].molecule_bonded_with;
					for(int j = 0; j < total_molecules; j++)
					{
						if((my_object_array[j].xaxis == x+1 && my_object_array[j].yaxis == y) && my_object_array[j].zaxis == z)
						{
							my_object_array[j].molecule_type = lattice[x][y][z].molecule_bonded_with;
                            break;	
						}
					}
                    lattice[x][y][z].molecule_bonded_with = 0;
                    lattice[x+1][y][z].molecule_bonded_with = 0;
                    lattice[x][y][z].bond_direction = 0;
                    lattice[x+1][y][z].bond_direction = 0;
                    for(int j = 0; j < total_molecules; j++)
					{
						if((my_object_array[j].xaxis == x && my_object_array[j].yaxis == y) && my_object_array[j].zaxis == z)
						{
                            my_object_array[j].molecule_bonded_with = 0;
                            my_object_array[j].bond_direction = 0;
                            break;	
						}
					}
                    for(int j = 0; j < total_molecules; j++)
					{
						if((my_object_array[j].xaxis == x+1 && my_object_array[j].yaxis == y) && my_object_array[j].zaxis == z)
						{
                            my_object_array[j].molecule_bonded_with = 0;
                            my_object_array[j].bond_direction = 0;
                            break;	
						}
					}
                    apoptosome_count --;
				}
				else if(lattice[x][y][z].bond_direction == 2)
				{				
					if(x > 0 && lattice[x-1][y][z] != null)
					{
						/*System.out.println("Lattice before");
						System.out.println("x y z");
						System.out.println(lattice[x][y][z].molecule_type);
						System.out.println(lattice[x][y][z].molecule_bonded_with);
						System.out.println(lattice[x][y][z].bond_direction);
						System.out.println("x-1 y z");
						System.out.println(lattice[x-1][y][z].molecule_type);
						System.out.println(lattice[x-1][y][z].molecule_bonded_with);
						System.out.println(lattice[x-1][y][z].bond_direction);
						System.out.println("MOA before");
						System.out.println("x y z");
						for(int j = 0; j < total_molecules; j++)
						{
							if(my_object_array[j].xaxis == x && my_object_array[j].yaxis == y && my_object_array[j].zaxis == z)
							{
								System.out.println(my_object_array[j].molecule_type);
	                            System.out.println(my_object_array[j].molecule_bonded_with);
	                            System.out.println(my_object_array[j].bond_direction);
	                            break;	
							}
						}
						System.out.println("x-1 y z");
						for(int j = 0; j < total_molecules; j++)
						{
							if(my_object_array[j].xaxis == x-1 && my_object_array[j].yaxis == y && my_object_array[j].zaxis == z)
							{
								System.out.println(my_object_array[j].molecule_type);
	                            System.out.println(my_object_array[j].molecule_bonded_with);
	                            System.out.println(my_object_array[j].bond_direction);
	                            break;	
							}
						}*/
						lattice[x][y][z].molecule_type = lattice[x-1][y][z].molecule_bonded_with;
						for(int j = 0; j < total_molecules; j++)
						{
							if(my_object_array[j].xaxis == x && my_object_array[j].yaxis == y && my_object_array[j].zaxis == z)
							{
								my_object_array[j].molecule_type = lattice[x-1][y][z].molecule_bonded_with;
	                            //my_object_array[j].molecule_bonded_with = 0;
	                            //my_object_array[j].bond_direction = 0;
	                            break;	
							}
						}
						lattice[x-1][y][z].molecule_type = lattice[x][y][z].molecule_bonded_with;
						for(int j = 0; j < total_molecules; j++)
						{
							if(my_object_array[j].xaxis == x-1 && my_object_array[j].yaxis == y && my_object_array[j].zaxis == z)
							{
								my_object_array[j].molecule_type = lattice[x][y][z].molecule_bonded_with;
	                            //my_object_array[j].molecule_bonded_with = 0;
	                            //my_object_array[j].bond_direction = 0;
	                            break;	
							}
						}
	                    lattice[x][y][z].molecule_bonded_with = 0;
	                    lattice[x-1][y][z].molecule_bonded_with = 0;
	                    lattice[x][y][z].bond_direction = 0;
	                    lattice[x-1][y][z].bond_direction = 0;
	                    apoptosome_count --;
						for(int j = 0; j < total_molecules; j++)
						{
							if(my_object_array[j].xaxis == x && my_object_array[j].yaxis == y && my_object_array[j].zaxis == z)
							{
								//my_object_array[j].molecule_type = lattice[x-1][y][z].molecule_bonded_with;
	                            my_object_array[j].molecule_bonded_with = 0;
	                            my_object_array[j].bond_direction = 0;
	                            break;	
							}
						}
						for(int j = 0; j < total_molecules; j++)
						{
							if(my_object_array[j].xaxis == x-1 && my_object_array[j].yaxis == y && my_object_array[j].zaxis == z)
							{
								//my_object_array[j].molecule_type = lattice[x][y][z].molecule_bonded_with;
	                            my_object_array[j].molecule_bonded_with = 0;
	                            my_object_array[j].bond_direction = 0;
	                            break;	
							}
						}
						/*System.out.println("Lattice after");
						System.out.println("x y z");
						System.out.println(lattice[x][y][z].molecule_type);
						System.out.println(lattice[x][y][z].molecule_bonded_with);
						System.out.println(lattice[x][y][z].bond_direction);
						System.out.println("x-1 y z");
						System.out.println(lattice[x-1][y][z].molecule_type);
						System.out.println(lattice[x-1][y][z].molecule_bonded_with);
						System.out.println(lattice[x-1][y][z].bond_direction);
						System.out.println("MOA after");
						System.out.println("x y z");
						for(int j = 0; j < total_molecules; j++)
						{
							if(my_object_array[j].xaxis == x && my_object_array[j].yaxis == y && my_object_array[j].zaxis == z)
							{
								System.out.println(my_object_array[j].molecule_type);
	                            System.out.println(my_object_array[j].molecule_bonded_with);
	                            System.out.println(my_object_array[j].bond_direction);
	                            break;	
							}
						}
						System.out.println("x-1 y z");
						for(int j = 0; j < total_molecules; j++)
						{
							if(my_object_array[j].xaxis == x-1 && my_object_array[j].yaxis == y && my_object_array[j].zaxis == z)
							{
								System.out.println(my_object_array[j].molecule_type);
	                            System.out.println(my_object_array[j].molecule_bonded_with);
	                            System.out.println(my_object_array[j].bond_direction);
	                            break;	
							}
						}*/
					}
				}
				else if(lattice[x][y][z].bond_direction == 3 && y != (lattice_dim-1) && lattice[x][y+1][z] != null)
				{
					lattice[x][y][z].molecule_type = lattice[x][y+1][z].molecule_bonded_with;
					for(int j = 0; j < total_molecules; j++)
					{
						if(my_object_array[j].xaxis == x && my_object_array[j].yaxis == y && my_object_array[j].zaxis == z)
						{
							my_object_array[j].molecule_type = lattice[x][y+1][z].molecule_bonded_with;
                            //my_object_array[j].molecule_bonded_with = 0;
                            //my_object_array[j].bond_direction = 0;
                            break;	
						}
					}
					lattice[x][y+1][z].molecule_type = lattice[x][y][z].molecule_bonded_with;
					for(int j = 0; j < total_molecules; j++)
					{
						if(my_object_array[j].xaxis == x && my_object_array[j].yaxis == y+1 && my_object_array[j].zaxis == z)
						{
							my_object_array[j].molecule_type = lattice[x][y][z].molecule_bonded_with;
                            //my_object_array[j].molecule_bonded_with = 0;
                            //my_object_array[j].bond_direction = 0;
                            break;	
						}
					}
                    lattice[x][y][z].molecule_bonded_with = 0;
                    lattice[x][y+1][z].molecule_bonded_with = 0;
                    lattice[x][y][z].bond_direction = 0;
                    lattice[x][y+1][z].bond_direction = 0;
                    apoptosome_count --;
					for(int j = 0; j < total_molecules; j++)
					{
						if(my_object_array[j].xaxis == x && my_object_array[j].yaxis == y && my_object_array[j].zaxis == z)
						{
							//my_object_array[j].molecule_type = lattice[x][y+1][z].molecule_bonded_with;
                            my_object_array[j].molecule_bonded_with = 0;
                            my_object_array[j].bond_direction = 0;
                            break;	
						}
					}
					for(int j = 0; j < total_molecules; j++)
					{
						if(my_object_array[j].xaxis == x && my_object_array[j].yaxis == y+1 && my_object_array[j].zaxis == z)
						{
							//my_object_array[j].molecule_type = lattice[x][y][z].molecule_bonded_with;
                            my_object_array[j].molecule_bonded_with = 0;
                            my_object_array[j].bond_direction = 0;
                            break;	
						}
					}
				}
				
				else if(lattice[x][y][z].bond_direction == 4)
				{
					if(y > 0 && lattice[x][y-1][z] != null)
					{
						lattice[x][y][z].molecule_type = lattice[x][y-1][z].molecule_bonded_with;
						for(int j = 0; j < total_molecules; j++)
						{
							if(my_object_array[j].xaxis == x && my_object_array[j].yaxis == y && my_object_array[j].zaxis == z)
							{
								my_object_array[j].molecule_type = lattice[x][y-1][z].molecule_bonded_with;
	                            //my_object_array[j].molecule_bonded_with = 0;
	                            //my_object_array[j].bond_direction = 0;
	                            break;	
							}
						}
						lattice[x][y-1][z].molecule_type = lattice[x][y][z].molecule_bonded_with;
						for(int j = 0; j < total_molecules; j++)
						{
							if(my_object_array[j].xaxis == x && my_object_array[j].yaxis == y-1 && my_object_array[j].zaxis == z)
							{
								my_object_array[j].molecule_type = lattice[x][y][z].molecule_bonded_with;
	                            //my_object_array[j].molecule_bonded_with = 0;
	                            //my_object_array[j].bond_direction = 0;
	                            break;	
							}
						}
	                    lattice[x][y][z].molecule_bonded_with = 0;
	                    lattice[x][y-1][z].molecule_bonded_with = 0;
	                    lattice[x][y][z].bond_direction = 0;
	                    lattice[x][y-1][z].bond_direction = 0;
	                    apoptosome_count --;
						for(int j = 0; j < total_molecules; j++)
						{
							if(my_object_array[j].xaxis == x && my_object_array[j].yaxis == y && my_object_array[j].zaxis == z)
							{
								//my_object_array[j].molecule_type = lattice[x][y-1][z].molecule_bonded_with;
	                            my_object_array[j].molecule_bonded_with = 0;
	                            my_object_array[j].bond_direction = 0;
	                            break;	
							}
						}
						for(int j = 0; j < total_molecules; j++)
						{
							if(my_object_array[j].xaxis == x && my_object_array[j].yaxis == y-1 && my_object_array[j].zaxis == z)
							{
								//my_object_array[j].molecule_type = lattice[x][y][z].molecule_bonded_with;
	                            my_object_array[j].molecule_bonded_with = 0;
	                            my_object_array[j].bond_direction = 0;
	                            break;	
							}
						}
					}
					
				}
				
				else if(lattice[x][y][z].bond_direction == 5 && z != (lattice_dim-1) && lattice[x][y][z+1] != null)
				{
					lattice[x][y][z].molecule_type = lattice[x][y][z+1].molecule_bonded_with;
					for(int j = 0; j < total_molecules; j++)
					{
						if(my_object_array[j].xaxis == x && my_object_array[j].yaxis == y && my_object_array[j].zaxis == z)
						{
							my_object_array[j].molecule_type = lattice[x][y][z+1].molecule_bonded_with;
                            //my_object_array[j].molecule_bonded_with = 0;
                            //my_object_array[j].bond_direction = 0;
                            break;	
						}
					}
					lattice[x][y][z+1].molecule_type = lattice[x][y][z].molecule_bonded_with;
					for(int j = 0; j < total_molecules; j++)
					{
						if(my_object_array[j].xaxis == x && my_object_array[j].yaxis == y && my_object_array[j].zaxis == z+1)
						{
							my_object_array[j].molecule_type = lattice[x][y][z].molecule_bonded_with;
                            //my_object_array[j].molecule_bonded_with = 0;
                            //my_object_array[j].bond_direction = 0;
                            break;	
						}
					}
                    lattice[x][y][z].molecule_bonded_with = 0;
                    lattice[x][y][z+1].molecule_bonded_with = 0;
                    lattice[x][y][z].bond_direction = 0;
                    lattice[x][y][z+1].bond_direction = 0;
                    apoptosome_count --;
					for(int j = 0; j < total_molecules; j++)
					{
						if(my_object_array[j].xaxis == x && my_object_array[j].yaxis == y && my_object_array[j].zaxis == z)
						{
							//my_object_array[j].molecule_type = lattice[x][y][z+1].molecule_bonded_with;
                            my_object_array[j].molecule_bonded_with = 0;
                            my_object_array[j].bond_direction = 0;
                            break;	
						}
					}
					for(int j = 0; j < total_molecules; j++)
					{
						if(my_object_array[j].xaxis == x && my_object_array[j].yaxis == y && my_object_array[j].zaxis == z+1)
						{
							//my_object_array[j].molecule_type = lattice[x][y][z].molecule_bonded_with;
                            my_object_array[j].molecule_bonded_with = 0;
                            my_object_array[j].bond_direction = 0;
                            break;	
						}
					}
				}
				else if(lattice[x][y][z].bond_direction == 6)
				{
					if(z > 0 && lattice[x][y][z-1] != null)
					{
						lattice[x][y][z].molecule_type = lattice[x][y][z-1].molecule_bonded_with;
						for(int j = 0; j < total_molecules; j++)
						{
							if(my_object_array[j].xaxis == x && my_object_array[j].yaxis == y && my_object_array[j].zaxis == z)
							{
								my_object_array[j].molecule_type = lattice[x][y][z-1].molecule_bonded_with;
	                            //my_object_array[j].molecule_bonded_with = 0;
	                            //my_object_array[j].bond_direction = 0;
	                            break;	
							}
						}
						lattice[x][y][z-1].molecule_type = lattice[x][y][z].molecule_bonded_with;
						for(int j = 0; j < total_molecules; j++)
						{
							if(my_object_array[j].xaxis == x && my_object_array[j].yaxis == y && my_object_array[j].zaxis == z-1)
							{
								my_object_array[j].molecule_type = lattice[x][y][z].molecule_bonded_with;
	                            //my_object_array[j].molecule_bonded_with = 0;
	                            //my_object_array[j].bond_direction = 0;
	                            break;	
							}
						}
	                    lattice[x][y][z].molecule_bonded_with = 0;
	                    lattice[x][y][z-1].molecule_bonded_with = 0;
	                    lattice[x][y][z].bond_direction = 0;
	                    lattice[x][y][z-1].bond_direction = 0;
	                    apoptosome_count --;
						for(int j = 0; j < total_molecules; j++)
						{
							if(my_object_array[j].xaxis == x && my_object_array[j].yaxis == y && my_object_array[j].zaxis == z)
							{
								//my_object_array[j].molecule_type = lattice[x][y][z-1].molecule_bonded_with;
	                            my_object_array[j].molecule_bonded_with = 0;
	                            my_object_array[j].bond_direction = 0;
	                            break;	
							}
						}
						for(int j = 0; j < total_molecules; j++)
						{
							if(my_object_array[j].xaxis == x && my_object_array[j].yaxis == y && my_object_array[j].zaxis == z-1)
							{
								//my_object_array[j].molecule_type = lattice[x][y][z].molecule_bonded_with;
	                            my_object_array[j].molecule_bonded_with = 0;
	                            my_object_array[j].bond_direction = 0;
	                            break;	
							}
						}
					}
				}
			}
		}
	}
	
	public static void main(String args[])
	{
		// By: Divyanshu Srivastava
		long start_time = System.nanoTime();
		Random r = new Random();
		double rand;
		int choice;
		for (int mc_runs=0; mc_runs<1; mc_runs++) {
			// This loop governs the total number of Monte_Carlo runs
			
			// Initialize the grid
			initialization();
			
			for (int mc_step=0; mc_step<mc_steps; mc_step++) {
				// This loop runs for each Monte_Carlo step.
				if (apoptosome_count == 2 && flag == 0)
					try {
						System.out.println("PAUSED");
						System.in.read();
						flag = 1;
						System.in.read();
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				if (dimer_count == 2 && flag2 == 0)
					try {
						System.out.println("PAUSED");
						System.in.read();
						flag2 = 1;
						System.in.read();
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				for (int mol=0; mol<total_molecules; mol++) {
					// This loops runs for each molecule in the lattice.
					rand = Math.random();
					choice = r.nextInt(total_molecules);
					if (rand < 0.5) {
						// Perform random walk in lattice
						if (my_object_array[choice].molecule_type == 1 || my_object_array[choice].molecule_type == 2)
							diffusion_cytoc_apaf(choice);
						else if (my_object_array[choice].molecule_type == 3) 
							if (r.nextDouble() < p_apoptosome_diffusion) diffusion_apoptosome(choice);
					} else {
						// Perform Reaction
						if (my_object_array[choice].molecule_bonded_with == 0) {
							// The randomly chosen molecule is not bonded
							bond_formation_by_cytoc_apaf(choice);
						} else {
							if (Math.random() < 0.5) {
								// Trying for bond breaking
								dissociation_apoptosome(choice);
							} else {
								// Trying for dimer formation
								bond_formation_by_apoptosome(choice);
							}
						}
					}
				}
				System.out.println(apoptosome_count + "\t" + dimer_count);
				System.out.println(mc_step);
				global_apoptosome_count[mc_step] = apoptosome_count;
				global_dimer_count[mc_step] = dimer_count;
			}			
			apoptosome_count = 0;
			dimer_count = 0;
			try
			{
			    PrintWriter pr = new PrintWriter("apoptosomes");    

			    for (int i=0; i<global_apoptosome_count.length ; i++)
			    {
			        pr.println(global_apoptosome_count[i]);
			    }
			    pr.close();
			}
			catch (Exception e)
			{
			    e.printStackTrace();
			    System.out.println("No such file exists.");
			}
			try
			{
			    PrintWriter pr = new PrintWriter("dimers");    

			    for (int i=0; i<global_dimer_count.length ; i++)
			    {
			        pr.println(global_dimer_count[i]);
			    }
			    pr.close();
			}
			catch (Exception e)
			{
			    e.printStackTrace();
			    System.out.println("No such file exists.");
			}
		}
		long stop_time = System.nanoTime();
		long execution_time = stop_time - start_time;
		System.out.println("Code Executed successfully");
		System.out.println("Execution Time: " + execution_time/1000000000 + " sec.");
	}
}