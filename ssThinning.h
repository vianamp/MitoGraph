#ifndef SSTHINNING_H
#define SSTHINNING_H

#include <cmath>
#include <cstdio>
#include <cstdlib>

	//===========================================================================
	//
	//   Thinning  algorithm  based  on the paper: "A 3D 6-subinteration thinning
	//   algorithm for extracting medial lines",  by  Kálman Palágyi  and  Attila
	//   Kuba.
	//
	//   Pattern Recognition Letters 19 (1998) 613-627.
	//   Kalman Palágyi email: palagyi@inf.u-szeged.hu.
	//
	// --------------------------------------------------------------------------
	//
	// Author: Matheus Palhares Viana
	// Email:       vianamp@gmail.com
	// Date:               06/04/2010
	//
	//===========================================================================

	/* 27-neighborhood considered in this routine:

						   x-1 x x+1

							0  1  2    z+1
	   y-1					3  4  5     z
							6  7  8    z-1
					 9 10 11
		y			12    13
					14 15 16
			17 18 19
	   y+1  20 21 22
			23 24 25

	*/

	//===========================================================================

    extern int ssdx[26], ssdy[26], ssdz[26];

	class ssMask {								 // Class  used  to  store a mask
												 // and its  rotations (90°, 180°
	  private:									 // and 270°).
		int *vector;
		int *vector90;
		int *vector180;
		int *vector270;
		char direction;
												 // Static functions to match and
												 // rotate the vectors

		static bool matchf(int ***Vol, int *vec);
		static void rotate(int vector[26],int *vector_rot,char axis);

	  public:
		ssMask();
		~ssMask();
		void print_mask();
		void print_masks();
		void gen_rotations();					 // Generates  the rotated  masks
		void print_direction();
		void set_direction(char d);
		bool match(int ***Vol);			 		 // Tests  if  one  of  the  four
												 // masks  (0°,90°,180° and 270°)
												 // matchs with the  neighborhood
												 // of the voxel p  in the volume
												 // Vol.

		void set_mask_from_u(int umask[26]);	 // Generates the four masks from
												 // the   mask  "umask"   in   up
												 // direction.
	};

	//===========================================================================

	class DirssMasks {							 // Class used to store the masks
												 // M1, M2, M3, M4, M5 and M6 for
												 // the direction "direction".
		private:								 // Observe that the i-th element
			ssMask *BaseMasks;					 // of   the   vector   BaseMasks
												 // corresponds  to  the  mask Mi
			public:								 // and this mask contains  their
			DirssMasks(char direction);			 // four  rotations.   Then,  the
			~DirssMasks();						 // class DirMask includes 6x4=24
			bool match(int ***Vol);		 		 // masks for the same direction.

	};

	//===========================================================================

	class ssThinVox {							 // This class is an interface to
												 // check  the  match for all six
												 // directions.
		private:
			DirssMasks *DirMask[6];

		public:
        
			ssThinVox();
			~ssThinVox();
			bool match(int direction, int ***V);

	};

#endif
