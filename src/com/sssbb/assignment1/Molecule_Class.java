package com.sssbb.assignment1;

public class Molecule_Class {
	
	int molecule_type;			// 1->Cytochrome_C 2->Apaf 3->Apoptosome 4->Dimer
	int molecule_bonded_with;	// 1->Cytochrome_C 2->Apaf
	int xaxis;					// x-coordinate in the lattice
	int yaxis;					// y-coordinate in the lattice
	int zaxis;					// z-coordinate in the lattice
	int bond_direction;			// This is intuitive
	
	public Molecule_Class(int molecule_type, int molecule_bonded_with, int xaxis, int yaxis, int zaxis, int bond_direction)
	{
		this.molecule_type = molecule_type;
		this.molecule_bonded_with = molecule_bonded_with;
		this.xaxis = xaxis;
		this.yaxis = yaxis;
		this.zaxis = zaxis;
		this.bond_direction = bond_direction;
	}
}
