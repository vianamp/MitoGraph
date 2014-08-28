	#include "ssThinning.h"


												 // Templates for 27-neighborhood
												 // and  the  six masks of the up
												 // direction  (paper:  fig.  6).

    int ssMu[6][26] = {{0,0,0,0,0,0,0,0,0,3,3,3,3,3,3,3,3,3,3,3,3,1,3,3,3,3},
                       {2,2,2,0,0,0,0,0,0,2,1,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2},
                       {2,2,2,0,0,2,0,0,2,2,1,2,2,1,2,2,2,2,2,2,2,1,2,2,2,2},
                       {0,0,1,0,0,0,0,0,0,2,2,1,2,2,2,2,2,2,2,2,2,1,2,2,2,2},
                       {0,0,0,0,0,0,0,0,0,3,3,3,3,3,0,0,0,3,1,3,3,0,3,0,0,0},
                       {0,0,0,0,0,0,0,0,0,2,2,2,0,2,0,0,2,2,1,2,0,0,1,0,0,2}};


	//===========================================================================

	ssMask::ssMask() {								 // Constructor
		   vector = new int[26];
		 vector90 = new int[26];
		vector180 = new int[26];
		vector270 = new int[26];
		direction = ' ';
		for (int i=26;i--;) vector[i] = vector90[i] = vector180[i] = vector270[i] = -1;
	}

	ssMask::~ssMask() {								 // Destructor
		delete[] vector;
		delete[] vector90;
		delete[] vector180;
		delete[] vector270;
	}

	void ssMask::set_direction(char d) { direction = d; }

	void ssMask::print_direction() { printf("%c\n",direction); }

	void ssMask::set_mask_from_u(int umask[26]) {// Uses  the umask in  order  to
		if (direction == ' ') return;            // generate  the masks for other
		rotate(umask,vector,'x');                // directions.     Indeed,    to
		if (direction == 'n') return;            // generate  the  "n" version of
		rotate(vector,vector,'x');               // mask M2, showed in he  figure
		if (direction == 'd') return;            // 3 of the paper, we just to do
		rotate(vector,vector,'x');               // a rotation  of M2  about  the
		if (direction == 's') return;            // axis  x.   To  perform   this
		rotate(vector,vector,'x');				 // transformation,     we   use:
		if (direction == 'u') return;	         // >>> rotate(umask,vector,'x'),
		rotate(vector,vector,'z');               // which rotates the mask  umask
		if (direction == 'w') return;            // about   "x"  and  stores  the
		rotate(vector,vector,'z');               // result in "vector".
		rotate(vector,vector,'z');
		if (direction == 'e') return;
	}

	void ssMask::gen_rotations() {				 // Rotations for the  directions
		char axis;								 // (u,d), (n,s) and  (w,e)  must
		if (direction == ' ') return;			 // be  performed  about the same
												 // axis: y, z and x, respec.

		if (direction == 'u' || direction == 'd') axis = 'y';
		if (direction == 'n' || direction == 's') axis = 'z';
		if (direction == 'w' || direction == 'e') axis = 'x';

			rotate(vector,vector90,axis);		 // 90° rotation
		 rotate(vector90,vector180,axis);		 // 90° rotation
		rotate(vector180,vector270,axis);		 // 90° rotation
	}

	void ssMask::print_mask() {
		for (int i=0;i<26;i++) printf("%d ",vector[i]);
		printf("\n");
	}

	void ssMask::print_masks() {
		print_direction();
		printf("0\t=\t");
		for (int i=0;i<26;i++) printf("%d ",vector[i]);
		printf("\n");
		printf("90\t=\t");
		for (int i=0;i<26;i++) printf("%d ",vector90[i]);
		printf("\n");
		printf("180\t=\t");
		for (int i=0;i<26;i++) printf("%d ",vector180[i]);
		printf("\n");
		printf("270\t=\t");
		for (int i=0;i<26;i++) printf("%d ",vector270[i]);
		printf("\n");
	}


	void ssMask::rotate(int vector[26],int *vector_rot,char axis) {
		int register i;
		int x, y, z, rt[9], aux[3][3][3];		 // Variables and auxiliar matrix
		int rx[9] = {1,0,0,0,0,1,0,-1,0};		 // Rotation matrices (Rx, Ry and
		int ry[9] = {0,0,-1,0,1,0,1,0,0};		 // Rz) are  stored  as  vectors.
		int rz[9] = {0,1,0,-1,0,0,0,0,1};		 // Note that Ri is the  rotation
												 // matrix about the axis i.
		for (i=9; i--;) {
			if(axis=='x') rt[i] = rx[i];
			if(axis=='y') rt[i] = ry[i];
			if(axis=='z') rt[i] = rz[i];
		}

		for (i=26; i--;) {						 // Performing  the  rotation  of
												 // the  vector   "vector".   The
												 // result   is   stored  in  the
												 // temporary matrix.

			x = rt[0]*ssdx[i] + rt[1]*ssdy[i] + rt[2]*ssdz[i];
			y = rt[3]*ssdx[i] + rt[4]*ssdy[i] + rt[5]*ssdz[i];
			z = rt[6]*ssdx[i] + rt[7]*ssdy[i] + rt[8]*ssdz[i];
			aux[x+1][y+1][z+1] = vector[i];
		}
												 // Updating  the result  vector.

		for (i=26; i--;) vector_rot[i] = aux[ssdx[i]+1][ssdy[i]+1][ssdz[i]+1];
	}

	// There are four types of indexes in the masks: 0, 1, 2 and 3. We need worry
	// about the indexes 0, 1 and 3 bacause 2 indicates "does not matter", ie the
	// value of this  voxel  in the  volume can be any.  Indexes 0 and 1  must be
	// exactly the same as in the  mask as in the volume. Index 3 is a particular
	// case. As described in the paper, at last one voxel of the volume with this
	// index must be 1.

	bool ssMask::matchf(int ***Vol, int *vec) {
		int v, q3 = 0;
		int register i;
		for (i=26;i--;) {
			v = Vol[1+ssdx[i]][1+ssdy[i]][1+ssdz[i]];
			if (vec[i]==0&&v!=0) return false;
			if (vec[i]==1&&v!=1) return false;
			if (vec[i]==3) q3 = (v==1||q3==2) ? 2 :1;
		}
												 // If  q3  remains  equal  to 0,
												 // there is no values 3  in  the
												 // the mask. If q3 == 2 there is
												 // value 3 in  the mask  and  at
												 // last one  voxel in the volume
												 // has value 1.
		if (q3==1) return false;
		return true;
	}

	bool ssMask::match(int ***Vol) {	 // Verify if the rotations match.
		   if (matchf(Vol,vector)) return true;
		 if (matchf(Vol,vector90)) return true;
		if (matchf(Vol,vector180)) return true;
		if (matchf(Vol,vector270)) return true;
		return false;
	}

	//===========================================================================

	DirssMasks::DirssMasks(char direction) {	 // Constructor: the 24 maks  are
		BaseMasks = new ssMask[6];				 // generated at this stage.
		for (int i=6;i--;) {
			BaseMasks[i].set_direction(direction);
			BaseMasks[i].set_mask_from_u(ssMu[i]);
			BaseMasks[i].gen_rotations();
		}
	}

	DirssMasks::~DirssMasks() { delete[] BaseMasks; }

	bool DirssMasks::match(int ***Vol) {
		for (int i=6;i--;) {
			if (BaseMasks[i].match(Vol)) return true;
		}
		return false;
	}

	//===========================================================================

	ssThinVox::ssThinVox() {					 // Constructor:
		DirMask[0] = new DirssMasks('u'); 		 // Generates the  24  masks  for
		DirMask[1] = new DirssMasks('d'); 		 // each  direction (6x24 = 144).
		DirMask[2] = new DirssMasks('n');
		DirMask[3] = new DirssMasks('s');
		DirMask[4] = new DirssMasks('e');
		DirMask[5] = new DirssMasks('w');
	}

	ssThinVox::~ssThinVox() {					 // Destructor
		delete DirMask[0];
		delete DirMask[1];
		delete DirMask[2];
		delete DirMask[3];
		delete DirMask[4];
		delete DirMask[5];
	}

	bool ssThinVox::match(int direction, int ***Vol) {
		if (DirMask[direction]->match(Vol)) return true;
		return false;
	}

	//===========================================================================
