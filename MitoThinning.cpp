// =====================================================================================================
// MitoThinning: This routine receives as input an ImageData binary from MitoGraph main routine a and
// apply a thinning process over it, resulting in a new but topologically equivalent image.
//
// Matheus P. Viana - vianamp@gmail.com - 2014.06.10
// Susanne Rafelski Lab, University of California Irvine
// =====================================================================================================

#include "ssThinning.h"
#include "MitoThinning.h"

    int ssdx[26] = {-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1};
    int ssdy[26] = {-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    int ssdz[26] = { 1, 1, 1, 0, 0, 0,-1,-1,-1, 1, 1, 1, 0, 0,-1,-1,-1, 1, 1, 1, 0, 0, 0,-1,-1,-1};

// Routine used to save .gnet and .coo files representing
// the skeleton of the mitochondrial network.
void ExportGraphFiles(vtkSmartPointer<vtkPolyData> PolyData, long int nnodes, const char Prefix[]);

// Routine to create the file _nodes.vtk. This file contains
// little speres located at the junctions (nodes) coordinates
// and might have the nodes label depending on if the variable
// export_node_labels is true or false.
void ExportNodes(vtkSmartPointer<vtkPolyData> PolyData, long int nnodes, long int *ValidId, const char Prefix[]);

// Routine used to check whether a voxel is locates at the border
// of the binary image or not. A border voxel is defined as those
// that have at least one of their 6-neighbors equal to zero.
bool IsBorder(vtkIdType id, vtkSmartPointer<vtkImageData> Image);

// Routine used to smooth the parametric curves that describe the
// skeleton edges. A simple average filter is used to do the job.
// Sigma represents the number of times the filter is applied.
void SmoothEdgesCoordinates(vtkSmartPointer<vtkPolyData> PolyData, double sigma);

// Estimate the mitochondrial volume by counting the number of
// pixels in the binary image used as input for thinning.
void GetVolumeFromVoxels(vtkSmartPointer<vtkImageData> Image, double *attributes);

// Estimate the mitochondrial volume by using the skeleton
// total length and assuming constant radius.
void GetVolumeFromSkeletonLength(vtkSmartPointer<vtkPolyData> PolyData, double *attributes);

// Calculate the length of a given edge.
double GetEdgeLength(vtkIdType edge, vtkSmartPointer<vtkPolyData> PolyData);

// Returns the number of voxels around the voxel (x,y,z) with
// value given different of "value".
char GetNumberOfNeighborsWithoutValue(vtkSmartPointer<vtkImageData> Image, int x, int y, int z, long int value);

// Returns the number of voxels around the voxel (x,y,z) with
// value different of "value" in the vector "Volume".
char GetNumberOfNeighborsWithoutValue(vtkSmartPointer<vtkImageData> Image, vtkSmartPointer<vtkLongArray> Volume, int x, int y, int z, long int value);

// Returns the number of voxels around the voxel (x,y,z) with
// value given by "value".
char GetNumberOfNeighborsWithValue(vtkSmartPointer<vtkImageData> Image, int x, int y, int z, long int value);

// Returns the number of voxels around the voxel (x,y,z) with
// value "value" in the vector "Volume".
char GetNumberOfNeighborsWithValue(vtkSmartPointer<vtkImageData> Image, vtkSmartPointer<vtkLongArray> Volume, int x, int y, int z, long int value);

// Returns one neighbor of (x,y,z) with value "value" in the
// vector "Volume".
vtkIdType GetOneNeighborWithValue(int x, int y, int z, vtkSmartPointer<vtkImageData> Image, vtkSmartPointer<vtkLongArray> Volume, long int value);

// Returns one neighbor of (x,y,z) with value different of
// "value" in the vector "Volume".
vtkIdType GetOneNeighborWithoutValue(int x, int y, int z, vtkSmartPointer<vtkImageData> Image, vtkSmartPointer<vtkLongArray> Volume, long int value);

// Returns one neighbor of (x,y,z) with value different of
// "value".
vtkIdType GetOneNeighborWithoutValue(int x, int y, int z, vtkSmartPointer<vtkImageData> Image, long int value);

// Merge together junctions that touch each other.
bool JunctionsMerge(std::list<vtkIdType> Junctions, vtkSmartPointer<vtkImageData> Image, vtkSmartPointer<vtkLongArray> Volume);

// Track an edge starting at voxel (x,y,z) in the volume "Volume".
std::list<vtkIdType> GetEdgeStartingAt(int x, int y, int z, vtkSmartPointer<vtkImageData> Image, vtkSmartPointer<vtkLongArray> Volume);

// Returns an adjacency edge "edge_label".
long int GetOneAdjacentEdge(vtkSmartPointer<vtkPolyData> PolyData, long int edge_label, long int junction_label, bool *found_on_left);

// Track all the nodes and edges of a 3D structured thinned by
// the routine Thinning3D.
vtkSmartPointer<vtkPolyData> Skeletonization(vtkSmartPointer<vtkImageData> Image, const char FileName[], double *attributes);

// Replace two edges A and B attached to a node of degree 2
// with one edge C that corresponds to A + B. The degree of the
// node is set to -1 and A and B are moreved from PolyData.
bool MergeEdgesOfDegree2Nodes(vtkSmartPointer<vtkPolyData> PolyData, int *K);

/* ================================================================
   I/O ROUTINES
=================================================================*/

void SaveImageData(vtkSmartPointer<vtkImageData> Image, const char FileName[], bool _resample) {
    #ifdef DEBUG
        printf("Saving ImageData File...\n");
    #endif

    vtkSmartPointer<vtkStructuredPointsWriter> writer = vtkSmartPointer<vtkStructuredPointsWriter>::New();

    if (_resample) {
        #ifdef DEBUG
            printf("\tResampling data...%f\t%f\n",_dxy,_dz);
        #endif
        vtkSmartPointer<vtkImageResample> Resample = vtkSmartPointer<vtkImageResample>::New();
        Resample -> SetInterpolationModeToLinear();
        Resample -> SetDimensionality(3);
        Resample -> SetInputData(Image);
        Resample -> SetAxisMagnificationFactor(0,1.0);
        Resample -> SetAxisMagnificationFactor(1,1.0);
        Resample -> SetAxisMagnificationFactor(2,_dz/_dxy);
        Resample -> Update();

        vtkSmartPointer<vtkImageData> ImageResampled = Resample -> GetOutput();
        ImageResampled -> SetSpacing(1,1,1);

        writer -> SetInputData(ImageResampled);
    } else {

        writer -> SetInputData(Image);

    }
    
    writer -> SetFileType(VTK_BINARY);
    writer -> SetFileName(FileName);
    writer -> Write();

    #ifdef DEBUG
        printf("\tFile Saved!\n");
    #endif
}

void SavePolyData(vtkSmartPointer<vtkPolyData> PolyData, const char FileName[], bool scale) {

    #ifdef DEBUG
        printf("Saving PolyData from XYZ list...\n");
    #endif

    #ifdef DEBUG
        printf("\t#Points in PolyData file: %llu.\n",(vtkIdType)PolyData->GetNumberOfPoints());
    #endif

    vtkSmartPointer<vtkPolyDataWriter> Writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    Writer -> SetFileType(VTK_BINARY);
    Writer -> SetFileName(FileName);
    Writer -> SetInputData(PolyData);
    Writer -> Write();

    #ifdef DEBUG
        printf("\tFile Saved!\n");
    #endif
}

void ExportGraphFiles(vtkSmartPointer<vtkPolyData> PolyData, long int nnodes, long int *ValidId, const char Prefix[]) {
    
    #ifdef DEBUG
        printf("Saving .coo file...\n");
    #endif

    char _fullpath[256];

    vtkPoints *Points = PolyData -> GetPoints();

    double r[3];
    long int node, exact_nnodes = 0;
    for (node = 0; node < nnodes; node++) {
        if (ValidId[node]>=0) exact_nnodes++;
    }

    sprintf(_fullpath,"%s.coo",Prefix);
    FILE *fcoo = fopen(_fullpath,"w");
    for (node = 0; node < nnodes; node++) {
        if (ValidId[node]>=0) {
            Points -> GetPoint(node,r);
            fprintf(fcoo,"%1.4f\t%1.4f\t%1.4f\n",_dxy*r[0],_dxy*r[1],_dz*r[2]);
        }
    }
    fclose(fcoo);

    double length;
    vtkIdType edge, npoints, i, j;
    sprintf(_fullpath,"%s.gnet",Prefix);
    FILE *fgnet = fopen(_fullpath,"w");
    fprintf(fgnet,"%ld\n",exact_nnodes);
    for (edge = 0; edge < PolyData -> GetNumberOfCells(); edge++) {
        npoints = PolyData -> GetCell(edge) -> GetNumberOfPoints();
        i = PolyData -> GetCell(edge) -> GetPointId(0);
        j = PolyData -> GetCell(edge) -> GetPointId(npoints-1);
        if ( ValidId[i] >= 0 && ValidId[j] >= 0 ) {
            length = GetEdgeLength(edge,PolyData);
            fprintf(fgnet,"%ld\t%ld\t%1.5f\n",ValidId[i],ValidId[j],length);
        }
    }
    fclose(fgnet);

}

void ExportNodes(vtkSmartPointer<vtkPolyData> PolyData, long int nnodes, long int *ValidId, const char Prefix[]) {
    
    #ifdef DEBUG
        if (_export_nodes_label) {
            printf("Exporting nodes and their labels...\n");
        } else {
            printf("Exporting nodes...\n");
        }
    #endif

    vtkPoints *Points = PolyData -> GetPoints();
    vtkSmartPointer<vtkAppendPolyData> Append = vtkSmartPointer<vtkAppendPolyData>::New();
    Append -> SetOutputPointsPrecision(vtkAlgorithm::DEFAULT_PRECISION);

    double r[3];
    long int node;
    char node_txt[16];
    for (node = 0; node < nnodes; node++) {
        if ( ValidId[node] >= 0 ) {
            Points -> GetPoint(node,r);
        
            vtkSmartPointer<vtkSphereSource> Node = vtkSmartPointer<vtkSphereSource>::New();
            if (_scale_polydata_before_save) {
                Node -> SetRadius(_dxy);
                Node -> SetCenter(_dxy*r[0],_dxy*r[1],_dz*r[2]);
            } else {
                Node -> SetRadius(2.0);
                Node -> SetCenter(r[0],r[1],r[2]);
            }

            Node -> SetThetaResolution(12);
            Node -> SetPhiResolution(12);
            Node -> Update();
            Append -> AddInputData(Node->GetOutput());
            Append -> Update();

            if (_export_nodes_label) {
                sprintf(node_txt,"%ld",ValidId[node]);
                vtkSmartPointer<vtkVectorText> PolyText = vtkSmartPointer<vtkVectorText>::New();
                PolyText -> SetText(node_txt);
                PolyText -> Update();

                vtkSmartPointer<vtkTransform> T = vtkSmartPointer<vtkTransform>::New();
                if (_scale_polydata_before_save) {
                    T -> Translate(_dxy*(r[0]+1),_dxy*(r[1]+1),_dz*r[2]);
                    T -> Scale(2*_dxy,2*_dxy,1);
                } else {
                    T -> Translate(r[0]+2,r[1],r[2]);
                    T -> Scale(2,2,1);
                }

                vtkSmartPointer<vtkTransformPolyDataFilter> Trans = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
                Trans -> SetInputData(PolyText->GetOutput());
                Trans -> SetTransform(T);
                Trans -> Update();

                Append -> AddInputData(Trans -> GetOutput());
                Append -> Update();
            }
        }

    }

    char _fullpath[256];
    sprintf(_fullpath,"%s_nodes.vtk",Prefix);
    
    SavePolyData(Append->GetOutput(),_fullpath,false);
}


/* ================================================================
   IMAGE TRANSFORMATION
=================================================================*/

void CleanImageBoundaries(vtkSmartPointer<vtkImageData> ImageData) {
    #ifdef DEBUG
        printf("Cleaning the image boundaries...\n");
    #endif
    int p, q;
    int *Dim = ImageData -> GetDimensions();
    for ( p = Dim[0]; p--; ) {
        for ( q = Dim[1]; q--; ) {
            ImageData -> SetScalarComponentFromDouble(p,q,0,0,0);
            ImageData -> SetScalarComponentFromDouble(p,q,Dim[2]-1,0,0);
        }
    }
    for ( p = Dim[0]; p--; ) {
        for ( q = Dim[2]; q--; ) {
            ImageData -> SetScalarComponentFromDouble(p,0,q,0,0);
            ImageData -> SetScalarComponentFromDouble(p,Dim[1]-1,q,0,0);
        }
    }
    for ( p = Dim[1]; p--; ) {
        for ( q = Dim[2]; q--; ) {
            ImageData -> SetScalarComponentFromDouble(0,p,q,0,0);
            ImageData -> SetScalarComponentFromDouble(Dim[0]-1,p,q,0,0);
        }
    }
}

/* ================================================================
   AUXILIAR ROUTINES
=================================================================*/


bool IsBorder(vtkIdType id, vtkSmartPointer<vtkImageData> Image) {
    double r[3];
    int x, y, z;
    Image -> GetPoint(id,r);
    x = (int)r[0]; y = (int)r[1]; z = (int)r[2];
    if(!Image->GetScalarComponentAsDouble(x,y,z,0))   return false;
    if(!Image->GetScalarComponentAsDouble(x+1,y,z,0)) return true;
    if(!Image->GetScalarComponentAsDouble(x-1,y,z,0)) return true;
    if(!Image->GetScalarComponentAsDouble(x,y+1,z,0)) return true;
    if(!Image->GetScalarComponentAsDouble(x,y-1,z,0)) return true;
    if(!Image->GetScalarComponentAsDouble(x,y,z+1,0)) return true;
    if(!Image->GetScalarComponentAsDouble(x,y,z-1,0)) return true;
    return false;
}

void SmoothEdgesCoordinates(vtkSmartPointer<vtkPolyData> PolyData, double sigma){

    double r1[3], r2[3], r3[3];
    long int id, line, nlines, n, s;
    vtkPoints *Points = PolyData -> GetPoints();

    vtkCell *Line;
    nlines = PolyData -> GetNumberOfCells();

    for (line = 0; line < nlines; line++){
        Line = PolyData -> GetCell(line);
        n = Line -> GetNumberOfPoints();
        double *X = new double[n];
        double *Y = new double[n];
        double *Z = new double[n];
        for (s = 0; s < int(sigma); s++) {
            for (id = 1; id < n-1; id++ ) {
                Points -> GetPoint(Line->GetPointId(id-1),r1);
                Points -> GetPoint(Line->GetPointId(id+0),r2);
                Points -> GetPoint(Line->GetPointId(id+1),r3);
                X[id] = 0.125*(r1[0]+6*r2[0]+r3[0]);
                Y[id] = 0.125*(r1[1]+6*r2[1]+r3[1]);
                Z[id] = (1.0/3.0)*(r1[2]+r2[2]+r3[2]);
            }
            for (id = 1; id < n-1; id++ ) {
                Points -> SetPoint(Line->GetPointId(id),X[id],Y[id],Z[id]);
            }
        }
        delete[] X; delete[] Y; delete[] Z;
    }
    PolyData -> Modified();
}

void GetVolumeFromVoxels(vtkSmartPointer<vtkImageData> Image, double *attributes) {
    double v;
    unsigned long int nv = 0;
    for (vtkIdType id=Image->GetNumberOfPoints();id--;) {
        v = Image -> GetPointData() -> GetScalars() -> GetTuple1(id);
        if (v) nv++;
    }
    attributes[0] = nv * (_dxy * _dxy * _dz);
}

void GetVolumeFromSkeletonLength(vtkSmartPointer<vtkPolyData> PolyData, double *attributes) {
    double r1[3], r2[3], length = 0.0;
    vtkPoints *Points = PolyData -> GetPoints();
    for (vtkIdType edge=PolyData->GetNumberOfCells();edge--;) {
        length += GetEdgeLength(edge,PolyData);
    }
    attributes[1] = length;
    attributes[2] = length * (acos(-1.0)*pow(_rad,2));
}

double GetEdgeLength(vtkIdType edge, vtkSmartPointer<vtkPolyData> PolyData) {
    double r1[3], r2[3];
    double length = 0.0;
    for (vtkIdType n = 1; n < PolyData->GetCell(edge)->GetNumberOfPoints(); n++) {
        PolyData -> GetPoint(PolyData->GetCell(edge)->GetPointId(n-1),r1);
        PolyData -> GetPoint(PolyData->GetCell(edge)->GetPointId(n  ),r2);
        length += sqrt(pow(_dxy*(r2[0]-r1[0]),2)+pow(_dxy*(r2[1]-r1[1]),2)+pow(_dz*(r2[2]-r1[2]),2));
    }
    return length;
}

/* ================================================================
   THINNING 3D
=================================================================*/

vtkSmartPointer<vtkPolyData> Thinning3D(vtkSmartPointer<vtkImageData> ImageData, const char FileName[], double *attributes) {

    #ifdef DEBUG
        char _imbinary[256];
        sprintf(_imbinary,"%s_binary.vtk",FileName);
        SaveImageData(ImageData,_imbinary);
    #endif

    vtkIdType N = ImageData -> GetNumberOfPoints();
    GetVolumeFromVoxels(ImageData,attributes);   //surface-based volume

    int x, y, z;
    double r[3], v, vl;
    vtkIdType ndels, id;
    ssThinVox *STV = new ssThinVox();

    CleanImageBoundaries(ImageData);

    int i, j, ***Vol = new int**[3];
    for (i = 0; i < 3; i++) {
        Vol[i] = new int*[3];
            for (j = 0; j < 3; j++) {
                Vol[i][j] = new int[3];
            }
    }

    #ifdef DEBUG
        printf("Starting thinning process...\n");
    #endif

    std::list<vtkIdType> ToBeDeleted;
    std::list<vtkIdType> OnTheSurface;
    std::list<vtkIdType>::iterator itId;

    do {
        ndels = 0;
        for (id = N; id--;) {
            if (IsBorder(id, ImageData)) {
                OnTheSurface.insert(OnTheSurface.begin(),id);
            }
        }

        for (int direction = 6; direction--;) {
            for ( itId=OnTheSurface.begin(); itId!=OnTheSurface.end(); itId++) {
                id = *itId;
                ImageData -> GetPoint(id,r);
                x = (int)r[0]; y = (int)r[1]; z = (int)r[2];
                v = ImageData -> GetScalarComponentAsDouble(x,y,z,0);
                if (v) {
                    for (i = 26; i--;) {
                        vl = ImageData -> GetScalarComponentAsDouble(x+ssdx[i],y+ssdy[i],z+ssdz[i],0);
                        Vol[1+ssdx[i]][1+ssdy[i]][1+ssdz[i]] = (vl) ? 1 : 0;
                    }
                    if ( STV -> match(direction,Vol) ) ToBeDeleted.insert(ToBeDeleted.begin(),id);
                }
            }
            itId = ToBeDeleted.begin();
            while (itId != ToBeDeleted.end()) {
                ndels++;
                ImageData -> GetPoint(*itId,r);
                x = (int)r[0]; y = (int)r[1]; z = (int)r[2];
                ImageData -> SetScalarComponentFromDouble(x,y,z,0,0);
                ToBeDeleted.erase(itId++);
            }
        }

        #ifdef DEBUG
            printf("\t#Surface = %llu / #Deletions = %llu\n",(long long int)OnTheSurface.size(),ndels);
        #endif

        OnTheSurface.clear();

    } while(ndels);
    
    delete STV;

    #ifdef DEBUG
        printf("Thinning done!\n");
    #endif

    OnTheSurface.clear();
    ToBeDeleted.clear();

    return Skeletonization(ImageData,FileName,attributes);

}

/* ================================================================
   SKELETONIZATION
=================================================================*/

char GetNumberOfNeighborsWithoutValue(vtkSmartPointer<vtkImageData> Image, int x, int y, int z, long int value) {
    double r[3];
    char nn = 0;
    for (char k = 26; k--;) {
        if (Image->GetScalarComponentAsDouble(x+ssdx[k],y+ssdy[k],z+ssdz[k],0) != value) nn++;
    }
    return nn;
}

char GetNumberOfNeighborsWithoutValue(vtkSmartPointer<vtkImageData> Image, vtkSmartPointer<vtkLongArray> Volume, int x, int y, int z, long int value) {
    double r[3];
    char nn = 0;
    vtkIdType idk;
    for (char k = 26; k--;) {
        idk = Image->FindPoint(x+ssdx[k],y+ssdy[k],z+ssdz[k]);
        if (Volume->GetTuple1(idk) != value) nn++;
    }
    return nn;
}

char GetNumberOfNeighborsWithValue(vtkSmartPointer<vtkImageData> Image, int x, int y, int z, long int value) {
    double r[3];
    char nn = 0;
    for (char k = 26; k--;) {
        if (Image->GetScalarComponentAsDouble(x+ssdx[k],y+ssdy[k],z+ssdz[k],0) == value) nn++;
    }
    return nn;
}

char GetNumberOfNeighborsWithValue(vtkSmartPointer<vtkImageData> Image, vtkSmartPointer<vtkLongArray> Volume, int x, int y, int z, long int value) {
    double r[3];
    char nn = 0;
    vtkIdType idk;
    for (char k = 26; k--;) {
        idk = Image -> FindPoint(x+ssdx[k],y+ssdy[k],z+ssdz[k]);
        if (Volume->GetTuple1(idk) == value) nn++;
    }
    return nn;
}

vtkIdType GetOneNeighborWithValue(int x, int y, int z, vtkSmartPointer<vtkImageData> Image, vtkSmartPointer<vtkLongArray> Volume, long int value) {
    vtkIdType idk;
    for (char k = 26; k--;) {
        idk = Image -> FindPoint(x+ssdx[k],y+ssdy[k],z+ssdz[k]);
        if (Volume -> GetTuple1(idk) == value) {
            return idk;
        }
    }
    return 0;  // We can do it because, by construction the voxel at id 0 should always be empty
}

vtkIdType GetOneNeighborWithoutValue(int x, int y, int z, vtkSmartPointer<vtkImageData> Image, vtkSmartPointer<vtkLongArray> Volume, long int value) {
    vtkIdType idk;
    for (char k = 26; k--;) {
        idk = Image -> FindPoint(x+ssdx[k],y+ssdy[k],z+ssdz[k]);
        if (Volume -> GetTuple1(idk) && Volume -> GetTuple1(idk) != value) {
            return idk;
        }
    }
    return 0;  // We can do it because, by construction the voxel at id 0 should always be empty
}

vtkIdType GetOneNeighborWithoutValue(int x, int y, int z, vtkSmartPointer<vtkImageData> Image, long int value) {
    vtkIdType idk;
    for (char k = 26; k--;) {
        idk = Image -> FindPoint(x+ssdx[k],y+ssdy[k],z+ssdz[k]);
        if (Image->GetScalarComponentAsDouble(x+ssdx[k],y+ssdy[k],z+ssdz[k],0)!=value) return idk;
    }
    return 0;  // We can do it because, by construction the voxel at id 0 should always be empty
}

bool JunctionsMerge(std::list<vtkIdType> Junctions, vtkSmartPointer<vtkImageData> Image, vtkSmartPointer<vtkLongArray> Volume) {
    char nk;
    double r[3];
    vtkIdType id;
    int x, y, z, k;
    bool _has_changed = false;
    std::list<vtkIdType>::iterator itId;
    long int junction_label, neigh_junction_label;
    for (itId=Junctions.begin(); itId!=Junctions.end(); itId++) {
        Image -> GetPoint(*itId,r);
        junction_label = Volume -> GetTuple1(*itId);
        x = (int)r[0]; y = (int)r[1]; z = (int)r[2];
        for (k = 0; k < 26; k++) {
            id = Image->FindPoint(x+ssdx[k],y+ssdy[k],z+ssdz[k]);
            neigh_junction_label = Volume -> GetTuple1(id);
            if (junction_label < neigh_junction_label) {
                _has_changed = true;
                Volume -> SetTuple1(id,junction_label);
            }
        }
    }
    Volume -> Modified();
    return _has_changed;
}

std::list<vtkIdType> GetEdgeStartingAt(int x, int y, int z, vtkSmartPointer<vtkImageData> Image, vtkSmartPointer<vtkLongArray> Volume) {
    double r[3];
    vtkIdType idk;
    std::list<vtkIdType> Edge;
    do {
        idk = GetOneNeighborWithValue(x,y,z,Image,Volume,-1);
        if (idk) {
            Volume -> SetTuple1(idk,0);
            Edge.insert(Edge.end(),idk);
            Image -> GetPoint(idk,r);
            x = (int)r[0]; y = (int)r[1]; z = (int)r[2];
        }
    } while(idk);
    return Edge;
}

long int GetOneAdjacentEdge(vtkSmartPointer<vtkPolyData> PolyData, long int edge_label, long int original_source, bool *_common_source) {
    vtkCell *Edge;
    long int source, target;
    for (long int edge = PolyData->GetNumberOfCells();edge--;) {
        if (edge != edge_label) {
            Edge = PolyData -> GetCell((vtkIdType)edge);
            source = Edge -> GetPointId(0);
            target = Edge -> GetPointId(Edge->GetNumberOfPoints()-1);
            if (original_source==source) {
                *_common_source = true;
                return edge;
            }
            if (original_source==target) {
                *_common_source = false;
                return edge;
            }
        }
    }
    return -1;
}

vtkSmartPointer<vtkPolyData> Skeletonization(vtkSmartPointer<vtkImageData> Image, const char FileName[], double *attributes) {

    #ifdef DEBUG
        char _imthinn[256];
        sprintf(_imthinn,"%s_thinned.vtk",FileName);
        SaveImageData(Image,_imthinn);
    #endif

    double r[3];
    int x, y, z;
    vtkIdType id;
    long int junction_label = 1;
    vtkIdType N = Image -> GetNumberOfPoints();
    std::list<vtkIdType> Junctions;
    std::list<vtkIdType>::iterator itId;

    // Inside this IF statement, isolated voxels and isolated pairs of voxels
    // are expanded. If this is not done, these voxels will not be detected.
    // This requires O(N).
    if (_improve_skeleton_quality) {

        #ifdef DEBUG
            printf("Improving skeletonization [step 1: isolated single and pair of voxels]\n");
        #endif

        char nn;
        double rn[3];
        vtkIdType idn;
        int xn, yn, zn;
        for (id = N; id--;) {
            Image -> GetPoint(id,r);
            x = (int)r[0]; y = (int)r[1]; z = (int)r[2];
            if (Image->GetScalarComponentAsDouble(x,y,z,0)) {
                nn = GetNumberOfNeighborsWithoutValue(Image,x,y,z,0);
                if (nn==0) {
                    // Expanding isolated voxel
                    Image -> SetScalarComponentFromDouble(x+1,y,z,0,255);
                    Image -> SetScalarComponentFromDouble(x-1,y,z,0,255);
                } else if (nn==1) {
                    idn = GetOneNeighborWithoutValue(x,y,z,Image,0);
                    Image -> GetPoint(idn,rn);
                    xn = (int)rn[0]; yn = (int)rn[1]; zn = (int)rn[2];
                    if (GetNumberOfNeighborsWithoutValue(Image,xn,yn,zn,0)==1) {
                        // Expanding isolated pair of voxels
                        Image -> SetScalarComponentFromDouble(x-(int)(rn[0]-r[0]),y-(int)(rn[1]-r[1]),z-(int)(rn[2]-r[2]),0,255);
                        Image -> SetScalarComponentFromDouble(xn+(int)(rn[0]-r[0]),yn+(int)(rn[1]-r[1]),zn+(int)(rn[2]-r[2]),0,255);
                    }                    
                }
            }
        }
    }

    vtkSmartPointer<vtkLongArray> Volume = vtkSmartPointer<vtkLongArray>::New();
    Volume -> SetNumberOfComponents(0);
    Volume -> SetNumberOfTuples(N);

    // Inside this IF statement, we label all connected components
    // and we search for those cc that don't have junctions. If any
    // is found, we force it to have a junction by adding its first
    // voxel to the list Junctions. This fixes the problem of not
    // detecting loop-shaped ccs.
    if (_improve_skeleton_quality) {
        #ifdef DEBUG
            printf("Improving skeletonization [step 2: verifying connected components]\n");
        #endif
        long int n_fixed = 0;
        std::vector<long int> CSz;
        long int cc, ncc = LabelConnectedComponents(Image,Volume,CSz,26,0);

        bool *HasJunction = new bool[ncc];
        for (cc = ncc; cc--;) HasJunction[cc] = 0;
        for (id = N; id--;) {
            Image -> GetPoint(id,r);
            x = (int)r[0]; y = (int)r[1]; z = (int)r[2];
            cc = (long int)Volume -> GetTuple1(id);
            if (cc) {
                if (GetNumberOfNeighborsWithValue(Image,Volume,x,y,z,cc) != 2) {
                    HasJunction[(unsigned long int)(-cc-1)] = true;
                }
            }
            Image -> SetScalarComponentFromDouble(x,y,z,0,-cc);
        }
        for (id = N; id--;) {
            Image -> GetPoint(id,r);
            x = (int)r[0]; y = (int)r[1]; z = (int)r[2];
            cc = (long int)Volume -> GetTuple1(id);
            if (cc) {
                if (!HasJunction[(unsigned long int)(-cc-1)]) {
                    HasJunction[(unsigned long int)(-cc-1)] = true;
                    Junctions.insert(Junctions.begin(),id);
                    Volume -> SetTuple1(id,junction_label);
                    junction_label++;
                    n_fixed++;
                } else {
                    Volume -> SetTuple1(id,-1);
                }
            } else {
                Volume -> SetTuple1(id,0);
            }
        }
        delete[] HasJunction;

        #ifdef DEBUG
            printf("\t#Components fixed = %ld\n",n_fixed);
        #endif

        Image -> GetPointData() -> GetScalars() -> Modified();

    } else {

        for (id = N; id--;) {
            if (Image -> GetPointData() -> GetScalars() -> GetTuple1(id)) {
                Volume -> SetTuple1(id,-1);
            } else {
                Volume -> SetTuple1(id,0);
            }
        }
        Volume -> Modified();

    }

    #ifdef DEBUG
        printf("Starting skeletonization...\n");
        printf("\tSearching for junctions...\n");
    #endif

    for (id = N; id--;) {
        Image -> GetPoint(id,r);
        x = (int)r[0]; y = (int)r[1]; z = (int)r[2];
        if (Image->GetScalarComponentAsDouble(x,y,z,0)) {
            if (GetNumberOfNeighborsWithoutValue(Image,x,y,z,0) != 2) {
                Junctions.insert(Junctions.begin(),id);
                Volume -> SetTuple1(id,junction_label);
                junction_label++;
            }
        }
    }
    Volume -> Modified();

    if (!Junctions.size()) {
        printf("Z-stack seems to be empty. Aborting...\n");
        return 0;
    }

    #ifdef DEBUG
        printf("\t#Junctions before merging = %ld\n",Junctions.size());
    #endif

    // MERGING: During the merging process, voxels belonging to
    // the same junction are merged together forming nodes.
    while (JunctionsMerge(Junctions,Image,Volume));

    // Listing all nodes label we have until this point
    std::list<long int> Labels;
    std::list<long int>::iterator itLabel;
    for (itId=Junctions.begin(); itId!=Junctions.end(); itId++) {
        Labels.insert(Labels.begin(),(long int)Volume->GetTuple1(*itId));
    }

    // UNIQUE of labels
    Labels.sort(); // must sort first
    Labels.unique();
    long int NumberOfNodes = Labels.size();

    #ifdef DEBUG
        printf("\t#Junctions (nodes) after merging = %ld\n",NumberOfNodes);
    #endif

    // COORDINATES of nodes
    // Vectors with static size to make the average calculation easy.
    long int node;
    double *X = new double[NumberOfNodes]; //Vector for x-coordinate
    double *Y = new double[NumberOfNodes]; //Vector for y-coordinate
    double *Z = new double[NumberOfNodes]; //Vector for z-coordinate
    double *S = new double[NumberOfNodes]; //Vector for junctions size
       int *K = new int[NumberOfNodes];    //Vector for junctions degree
    for (node = 0; node < NumberOfNodes; node++) {
        K[node] = 0;
        X[node] = Y[node] = Z[node] = S[node] = 0.0;
    }

    for (itId=Junctions.begin(); itId!=Junctions.end(); itId++) {
        junction_label = (long int)Volume -> GetTuple1(*itId);
        itLabel = std::find(Labels.begin(),Labels.end(),junction_label);
        node = (long int)std::distance(Labels.begin(),itLabel);
        Volume -> SetTuple1(*itId,node+1);
        Image -> GetPoint(*itId,r);
        X[node] += r[0];
        Y[node] += r[1];
        Z[node] += r[2];
        S[node] ++;
    }
    Volume -> Modified();

    // Dynamic list to make easy the insertion of the coordinates
    // of remaining voxels detected during the edge tracking process.
    std::list<double> PointsListX;
    std::list<double> PointsListY;
    std::list<double> PointsListZ;
    for (node=0;node<NumberOfNodes;node++) {
        PointsListX.insert(PointsListX.end(),X[node]/S[node]);
        PointsListY.insert(PointsListY.end(),Y[node]/S[node]);
        PointsListZ.insert(PointsListZ.end(),Z[node]/S[node]);
    }

    delete[] X; delete[] Y; delete[] Z; delete[] S;

    // CellArray to add new edges as they are being tracked.
    vtkSmartPointer<vtkCellArray> EdgeArray = vtkSmartPointer<vtkCellArray>::New();

    bool _should_add;
    itId = Junctions.begin();
    long int voxel_label = NumberOfNodes;
    long int junction_label_left, junction_label_right;
    long int source_node, target_node;
    while (itId != Junctions.end()) {
        Image -> GetPoint(*itId,r);
        x = (int)r[0]; y = (int)r[1]; z = (int)r[2];
        if (GetNumberOfNeighborsWithValue(Image,Volume,x,y,z,-1)) {
            _should_add = true;

            //Tracking new edge
            std::list<vtkIdType> Edge = GetEdgeStartingAt(x,y,z,Image,Volume);

            //Identifying junctions on left side of the edge
            source_node = Volume -> GetTuple1(*itId);

            //Identifying junctions on right side of the edge
            Image -> GetPoint(Edge.back(),r);
            x = (int)r[0]; y = (int)r[1]; z = (int)r[2];
            target_node = Volume -> GetTuple1(GetOneNeighborWithoutValue(x,y,z,Image,Volume,source_node));

            // When target_node is equal to 0 at this point, we are
            // dealing with loops. These loops can be real loops or noise
            // loops of very small length. This very short loops are often
            // related to noise around junctions or 1-voxel holes in the
            // data.
            // @@PARAMETER: Minimum loop length (in voxels)
            if (!target_node) {
                target_node = source_node;
                if (Edge.size() < 3) _should_add = false;
            }

            if (_should_add) {
                // Creating new cell for this edge
                EdgeArray -> InsertNextCell(Edge.size()+2);
                // Adding node 1
                EdgeArray -> InsertCellPoint(source_node-1);
                // Addint the new points coordinates as well as the edge itself
                for (std::list<vtkIdType>::iterator itIde=Edge.begin(); itIde!=Edge.end(); itIde++) {
                    Image -> GetPoint(*itIde,r);
                    PointsListX.insert(PointsListX.end(),r[0]);
                    PointsListY.insert(PointsListY.end(),r[1]);
                    PointsListZ.insert(PointsListZ.end(),r[2]);
                    EdgeArray -> InsertCellPoint(voxel_label);
                    voxel_label++;
                }
                // Adding node 2
                EdgeArray -> InsertCellPoint(target_node-1);
                // Updating nodes degree
                K[source_node-1]++;
                K[target_node-1]++;
            }
            
        } else {
            Junctions.erase(itId++);
        }
    }

    vtkSmartPointer<vtkPoints> Points = vtkSmartPointer<vtkPoints>::New();
    Points -> SetNumberOfPoints(PointsListX.size());

    #ifdef DEBUG
        printf("#Points in vtkPoints = %lld\n",Points -> GetNumberOfPoints());
    #endif

    voxel_label = 0;
    std::list<double>::iterator itX = PointsListX.begin();
    std::list<double>::iterator itY = PointsListY.begin();
    std::list<double>::iterator itZ = PointsListZ.begin();
    while (itX != PointsListX.end()) {
        Points -> SetPoint(voxel_label,*itX,*itY,*itZ);
        itX++; itY++; itZ++;
        voxel_label++;
    }
    Points -> Modified();

    PointsListX.clear();
    PointsListY.clear();
    PointsListZ.clear();

    // Creating raw polyData
    vtkSmartPointer<vtkPolyData> PolyData = vtkPolyData::New();
    PolyData -> SetPoints(Points);
    PolyData -> SetLines(EdgeArray);
    PolyData -> Modified();
    PolyData -> BuildLinks();

    #ifdef DEBUG
        printf("\t#Edges before filtering = %lld\n",PolyData->GetNumberOfCells());
    #endif

    // PolyData filtering by removing degree-2 nodes. These nodes rise
    // for two reasons: 1) very short edges that are not detected,
    // although the bifurcation is detected. 2) Short loops that were
    // removed.
    while (MergeEdgesOfDegree2Nodes(PolyData,K));

    //Creating a list with final Ids of nodes after degree-2 nodes
    //deletetion.
    long int valid_id = 0;
    long int *ValidId = new long int[NumberOfNodes];
    for (node = 0; node < NumberOfNodes; node++) {
        if (K[node] > 0) {
            ValidId[node] = valid_id;
            valid_id++;
        } else {
            ValidId[node] = -1;
        }
    }

    #ifdef DEBUG
        printf("\t#Edges after filtering = %lld\n",PolyData->GetNumberOfCells());
        printf("Skeletonization done!\n");
    #endif

    SmoothEdgesCoordinates(PolyData,3);

    if (_export_graph_files) ExportGraphFiles(PolyData,NumberOfNodes,ValidId,FileName);

    ExportNodes(PolyData,NumberOfNodes,ValidId,FileName);

    GetVolumeFromSkeletonLength(PolyData,attributes); // total length and skeleton-length volume

    return PolyData;
}

bool MergeEdgesOfDegree2Nodes(vtkSmartPointer<vtkPolyData> PolyData, int *K) {

    // Given this edge: (source) o---->----o (target), it's necessary
    // to check whether either source or target are nodes of degree 2.

    std::list<vtkIdType> Merg;
    std::list<vtkIdType>::iterator itId;

    vtkCell *Edge, *NeighEdge;
    long int i, source, target, edge, edge_length, neigh_edge, neigh_edge_length;
    bool _common_source, _merged = false;

    for (edge = 0; edge < PolyData -> GetNumberOfCells(); edge++) {

        Edge = PolyData -> GetCell(edge);
        edge_length = Edge -> GetNumberOfPoints();
        source = (long int)Edge -> GetPointId(0);
        target = (long int)Edge -> GetPointId(Edge->GetNumberOfPoints()-1);

        if (K[source]==2 && source!=target) {

            //  [  neigh  | original ]
            //  o----?----o---->-----o
            //            s          t

            _merged = true;
            neigh_edge = GetOneAdjacentEdge(PolyData,edge,source,&_common_source);
            NeighEdge = PolyData -> GetCell((vtkIdType)neigh_edge);
            neigh_edge_length = NeighEdge -> GetNumberOfPoints();

            if (_common_source) {//  o----<----o---->----o
                for (i=0;i<edge_length;i++)
                    Merg.insert(Merg.begin(),PolyData->GetCell(edge)->GetPointId(i));
                for (i=1;i<neigh_edge_length;i++)
                    Merg.insert(Merg.end(),PolyData->GetCell(neigh_edge)->GetPointId(i));
            } else { //  o---->----o---->----o
                for (i=0;i<neigh_edge_length;i++)
                    Merg.insert(Merg.end(),PolyData->GetCell(neigh_edge)->GetPointId(i));                
                for (i=1;i<edge_length;i++)
                    Merg.insert(Merg.end(),PolyData->GetCell(edge)->GetPointId(i));
            }
            K[source] = -1; // Tagged as not valid node

        } else if (K[target]==2 && source!=target) {

            //  [ original |  neigh  ]
            //  o---->-----o----?----o
            //  s          t

            _merged = true;
            neigh_edge = GetOneAdjacentEdge(PolyData,edge,target,&_common_source);
            NeighEdge = PolyData -> GetCell((vtkIdType)neigh_edge);
            neigh_edge_length = NeighEdge -> GetNumberOfPoints();

            if (_common_source) {//  o---->----o---->----o
                for (i=0;i<edge_length;i++)
                    Merg.insert(Merg.end(),PolyData->GetCell(edge)->GetPointId(i));
                for (i=1;i<neigh_edge_length;i++)
                    Merg.insert(Merg.end(),PolyData->GetCell(neigh_edge)->GetPointId(i));
            } else {//  o---->----o----<----o
                for (i=0;i<edge_length;i++)
                    Merg.insert(Merg.end(),PolyData->GetCell(edge)->GetPointId(i));
                for (i=neigh_edge_length-1;i--;)
                    Merg.insert(Merg.end(),PolyData->GetCell(neigh_edge)->GetPointId(i));                
            }
            K[target] = -1; // Tagged as not valid node
        }

        if (_merged) {
            i = 0;
            // Removing the two pieces that were merged together
            PolyData -> BuildLinks();
            PolyData -> DeleteCell(edge);
            PolyData -> DeleteCell(neigh_edge);
            PolyData -> RemoveDeletedCells();
            // Adding the new "merged" edge
            vtkIdType IdList[Merg.size()];
            for (itId=Merg.begin(); itId!=Merg.end(); itId++) {
                IdList[i] = *itId;
                i++;
            }
            PolyData -> InsertNextCell(VTK_POLY_LINE,Merg.size(),IdList);
            Merg.clear();
            return true;
        }

    }
    return false;
}

long int LabelConnectedComponents(vtkSmartPointer<vtkImageData> ImageData, vtkSmartPointer<vtkDataArray> Volume, std::vector<long int> &CSz, int ngbh, double threshold) {

    #ifdef DEBUG
        printf("\tCalculating connected components...\n");
    #endif

    int *Dim = ImageData -> GetDimensions();

    vtkIdType i, s, ido, id;

    int x, y, z;
    double v, r[3];
    bool find = true;
    long long int ro[3];
    long int scluster, label;
    ro[0] = Dim[0] * Dim[1] * Dim[2];
    ro[1] = Dim[0] * Dim[1];

    vtkSmartPointer<vtkIdList> CurrA = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> NextA = vtkSmartPointer<vtkIdList>::New();

    Volume -> CopyComponent(0,ImageData->GetPointData()->GetScalars(),0);
    Volume -> Modified();

    label = 0;
    while (find) {
        for (s = 0; s < CurrA->GetNumberOfIds(); s++) {
            ido = CurrA -> GetId(s);
            ImageData -> GetPoint(ido,r);
            x = (int)r[0]; y = (int)r[1]; z = (int)r[2];
            for (i = 0; i < ngbh; i++) {
                id = ImageData -> FindPoint(x+ssdx_sort[i],y+ssdy_sort[i],z+ssdz_sort[i]);
                v = Volume -> GetTuple1(id);
                if (v > threshold) {
                    NextA -> InsertNextId(id);
                    Volume -> SetTuple1(id,-label);
                    scluster++;
                }
            }
        }
        if (!NextA->GetNumberOfIds()) {
            find = false;
            for (id=ro[0]; id--;) {
                v = Volume -> GetTuple1(id);
                if (v > threshold) {
                    find = true;
                    ro[0] = id;
                    break;
                }
            }
            if (label) {
                CSz.push_back(scluster);
            }
            if (find) {
                label++;
                scluster = 1;
                Volume -> SetTuple1(id,-label);
                CurrA -> InsertNextId(id);
            }
        } else {
            CurrA -> Reset();
            CurrA -> DeepCopy(NextA);
            NextA -> Reset();
        }
    }

    #ifdef DEBUG
        printf("\tNumber of detected components: %ld\n",(long int)CSz.size());
    #endif

    //return (long int)CSz->GetNumberOfTuples();
    return (long int)CSz.size();
}
