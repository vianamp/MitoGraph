// ==============================================================
// MitoGraph: Quantifying Mitochondrial Content in Live Cells
// Developed by Matheus P. Viana - vianamp@gmail.com - 2014.05.28
// Susanne Rafelski Lab, University of California Irvine
// Please, check the documentation at ?
// ==============================================================

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vtkMath.h>
#include <vtkImageFFT.h>
#include <vtkImageData.h>
#include <vtkDataArray.h>
#include <vtkDataObject.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredPoints.h>
#include <vtkInformationVector.h>
#include <vtkUnsignedCharArray.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkImageExtractComponents.h>
#include <vtkStructuredPointsReader.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkPolyDataWriter.h>
#include <vtkContourFilter.h>
#include <vtkDoubleArray.h>
#include <vtkTIFFReader.h>
#include <vtkPointData.h>
#include <vtkImageRFFT.h>
#include <vtkImageCast.h>

//#define DEBUG

// This routine returns the x of the id-th point of a 3D volume
// of size Dim[0]xDim[1]xDim[2]
int  GetX(unsigned long int id, int *Dim);

// This routine returns the y of the id-th point of a 3D volume
// of size Dim[0]xDim[1]xDim[2]
int  GetY(unsigned long int id, int *Dim);

// This routine returns the z of the id-th point of a 3D volume
// of size Dim[0]xDim[1]xDim[2]
int  GetZ(unsigned long int id, int *Dim);

// This routine returns the id of a point located at coordinate
// (x,y,z) of a 3D volume of size Dim[0]xDim[1]xDim[2]
unsigned long int GetId(int x, int y, int z, int *Dim);

// Get the reflected id of a point located at coordinate
// (x,y,z) of a 3D volume of size Dim[0]xDim[1]xDim[2]. The
// reflection is made with respect to the center of the
// 3D volume. This routine is used to shift the Gaussian
// derivatives to the corner of the volume before applying
// the Fourier transform
unsigned long int GetReflectedId(int x, int y, int z, int *Dim);

// Swap values
void Swap(double *x, double *y);

// Simple sorting algorithm and the output is such that
// l3 >= l2 >= l3
void Sort(double *l1, double *l2, double *l3);

// Calculate the Frobenius norm of a given 3x3 matrix
// http://mathworld.wolfram.com/FrobeniusNorm.html
double FrobeniusNorm(double M[3][3]);

// This routine fills a 3D volume (Kernel) with the a derivative
// of the gaussian function with variance sigma
// Use:
//      derivative = 1 for d2G/d2x
//      derivative = 2 for d2G/d2y
//      derivative = 3 for d2G/d2z
//      derivative = 4 for d2G/dxdy
//      derivative = 5 for d2G/dxdz
//      derivative = 6 for d2G/dydz
//      derivative = 0 for G (fills with the Gaussian function)
void GetGaussianDerivativeKernel(int derivative, vtkImageData *Kernel, double sigma);

// This routine uses the concept of scale space to calculate the
// derivatives of a given 3D volume. The Fourier transform is
// used to perform the convolution of the original image with
// the Gaussian derivatives
// ==============================================================
// @@FIXME: Try to use the complex multiplication from VTK
// ==============================================================
void GetImageDerivativeFourier(int derivative, double sigma, vtkImageData *ImageData, vtkDoubleArray *D);

// This routine calculate the Hessian matrix for each point
// of a 3D volume and its eigenvalues (Fourier Tranform Approach)
// ==============================================================
// @@FIXME: Ths function is producing negative values of vesselness
// ==============================================================
void GetHessianEigenvaluesFourier(double sigma, vtkImageData *Image, vtkDoubleArray *L1, vtkDoubleArray *L2, vtkDoubleArray *L3);

// Calculate the vesselness at each point of a 3D volume based
// based on the Hessian eigenvalues
void GetVesselness(double sigma, vtkImageData *Image, vtkDoubleArray *L1, vtkDoubleArray *L2, vtkDoubleArray *L3);

// This routine calculates the divergent filter of a 3D volume
// based on the orientation of the gradient vector field
void GetDivergentFilter(int *Dim, vtkDoubleArray *Scalars);

// This routine saves a 3D volume as VTK legacy file
void SaveImageData(vtkImageData *Image);

// This routine converts 16-bit volumes into 8-bit volumes by
// linearly scaling the original range of intensities [min,max]
// in [0,255] (http://rsbweb.nih.gov/ij/docs/guide/146-28.html)
vtkImageData *Convert16To8bit(vtkImageData *Image16);

// This routine saves a 3D polydata as VTK legacy file
void SavePolyData(vtkPolyData *PolyData, const char FileName[]);

// Given a 3D volume of size Dim[0]xDim[1]xDim[2], this routine
// adds a periodic border of size 0.1Dim[i] to the i-th dimension
// ==============================================================
// @@FIXME: Note that the corners of the expanded volume are not
// being filled
// ==============================================================
vtkImageData *AddBorder(vtkImageData *ImageData);

// This routine removes the border added by the routine above.
vtkImageData *RemoveBorder(int *DimB, vtkDoubleArray *BScalars, int *DimO);

// This routine uses a discrete differential operator to
// calculate the derivatives of a given 3D volume
void GetImageDerivativeDiscrete(vtkDataArray *Image, int *dim, char direction, vtkFloatArray *Derivative);

// This routine calculate the Hessian matrix for each point
// of a 3D volume and its eigenvalues (Discrete Approach)
void GetHessianEigenvaluesDiscrete(double sigma, vtkImageData *Image, vtkDoubleArray *L1, vtkDoubleArray *L2, vtkDoubleArray *L3);

/**========================================================
 Auxiliar functions
 =========================================================*/

int GetX(unsigned long int id, int *Dim) {
    return (int) id%Dim[0];
}

int GetY(unsigned long int id, int *Dim) {
    return (int)(id%(Dim[0]*Dim[1]))/Dim[0];
}

int GetZ(unsigned long int id, int *Dim) {
    return (int) id/(Dim[0]*Dim[1]);
}

unsigned long int GetId(int x, int y, int z, int *Dim) {
    return x+y*Dim[0]+z*Dim[0]*Dim[1];
}

unsigned long int GetReflectedId(int x, int y, int z, int *Dim) {
    int rx = ceil(0.5*Dim[0]);
    int ry = ceil(0.5*Dim[1]);
    int rz = ceil(0.5*Dim[2]);
    int sx = (x-(rx-0.5)<0) ? -rx : Dim[0]-rx;
    int sy = (y-(ry-0.5)<0) ? -ry : Dim[1]-ry;
    int sz = (z-(rz-0.5)<0) ? -rz : Dim[2]-rz;
    return GetId(x-sx,y-sy,z-sz,Dim);
}

void Swap(double *x, double *y) {
    double t = *y; *y = *x; *x = t;
}
void Sort(double *l1, double *l2, double *l3) {
    if (fabs(*l1) > fabs(*l2)) Swap(l1,l2);
    if (fabs(*l2) > fabs(*l3)) Swap(l2,l3);
    if (fabs(*l1) > fabs(*l2)) Swap(l1,l2);
}

double FrobeniusNorm(double M[3][3]) {
    double f = 0.0;
    for (int i = 3;i--;)
        for (int j = 3;j--;)
            f += M[i][j]*M[i][j];
    return sqrt(f);
}

/* ================================================================
   I/O ROUTINES
=================================================================*/

void SaveImageData(vtkImageData *Image, const char FileName[]) {
    #ifdef DEBUG
        printf("Saving ImageData File...\n");
    #endif

    vtkSmartPointer<vtkStructuredPointsWriter> writer = vtkSmartPointer<vtkStructuredPointsWriter>::New();
    #if (VTK_MAJOR_VERSION==5)
        writer->SetInput(Image);
    #else
        writer->SetInputData(Image);
    #endif
    writer->SetFileName(FileName);
    writer->Write();

    #ifdef DEBUG
        printf("File Saved!\n");
    #endif
}

void SavePolyData(vtkPolyData *PolyData, const char FileName[]) {
    #ifdef DEBUG
        printf("Saving PolyData File...\n");
    #endif

    vtkSmartPointer<vtkPolyDataWriter> PolyDataWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
    #if (VTK_MAJOR_VERSION==5)    
        PolyDataWriter -> SetInput(PolyData);
    #else
        PolyDataWriter -> SetInputData(PolyData);
    #endif
    PolyDataWriter -> SetFileName(FileName);
    PolyDataWriter -> Write();

    #ifdef DEBUG
        printf("File Saved!\n");
    #endif
}

/* ================================================================
   IMAGE TRANSFORM
=================================================================*/

vtkImageData *Convert16To8bit(vtkImageData *Image) {

    // 8-Bit images
    if (Image -> GetScalarType() == VTK_UNSIGNED_CHAR) {

        return Image;

    // 16-Bit images
    } else if (Image -> GetScalarType() == VTK_UNSIGNED_SHORT) {

        printf("Converting from 16-bit to 8-bit...\n");

        vtkImageData *Image8 = vtkImageData::New();
        Image8 -> ShallowCopy(Image);

        vtkDataArray *ScalarsShort = Image -> GetPointData() -> GetScalars();
        unsigned long int N = ScalarsShort -> GetNumberOfTuples();
        double range[2];
        ScalarsShort -> GetRange(range);

        printf("Original intensities range: [%d-%d]\n",(int)range[0],(int)range[1]);

        vtkSmartPointer<vtkUnsignedCharArray> ScalarsChar = vtkSmartPointer<vtkUnsignedCharArray>::New();
        ScalarsChar -> SetNumberOfComponents(1);
        ScalarsChar -> SetNumberOfTuples(N);
        
        double x, y;
        unsigned long int register id;
        for ( id = N; id--; ) {
            x = ScalarsShort -> GetTuple1(id);
            y = 255.0 * (x-range[0]) / (range[1]-range[0]);
            ScalarsChar -> SetTuple1(id,(unsigned char)y);
        }
        ScalarsChar -> Modified();

        Image8 -> GetPointData() -> SetScalars(ScalarsChar);
        return Image8;

    // Other depth
    } else {
        return NULL;
    }
}


/* ================================================================
   BORDER ROUTINES
=================================================================*/

vtkImageData *AddBorder(vtkImageData *ImageData) {
    #ifdef DEBUG
        printf("Adding border...\n");
    #endif

    int *NewDim = new int[3];
    int *Dim = ImageData -> GetDimensions();
    int N = ImageData -> GetNumberOfPoints();
    int dx = int(0.1*Dim[0]);
    int dy = int(0.1*Dim[1]);
    int dz = int(0.1*Dim[2]);
    int Lx = Dim[0] + 2*dx; Lx = (Lx%2) ? Lx+1 : Lx;
    int Ly = Dim[1] + 2*dy; Ly = (Ly%2) ? Ly+1 : Ly;
    int Lz = Dim[2] + 2*dz; Lz = (Lz%2) ? Lz+1 : Lz;
    NewDim[0] = Lx; NewDim[1] = Ly;  NewDim[2] = Lz;

    vtkDataArray *OScalars = ImageData -> GetPointData() -> GetScalars();
    vtkDoubleArray *BScalars = vtkDoubleArray::New();
    BScalars -> SetNumberOfTuples(Lx*Ly*Lz);
    BScalars -> FillComponent(0,0);

    double v;
    int x, y, z;
    unsigned long int id;
    for (x=0;x<dx;x++) {
        for (y=0;y<Dim[1];y++) {
            for (z=0;z<Dim[2];z++) {
                v = OScalars -> GetTuple1(GetId(x,y,z,Dim));
                BScalars -> SetTuple1(GetId(x+dx+Dim[0],y+dy,z+dz,NewDim),v);
                v = OScalars -> GetTuple1(GetId(Dim[0]-dx+x,y,z,Dim));
                BScalars -> SetTuple1(GetId(x,y+dy,z+dz,NewDim),v);
            }
        }
    }
    for (x=0;x<Dim[0];x++) {
        for (y=0;y<dy;y++) {
            for (z=0;z<Dim[2];z++) {
                v = OScalars -> GetTuple1(GetId(x,y,z,Dim));
                BScalars -> SetTuple1(GetId(x+dx,y+dy+Dim[1],z+dz,NewDim),v);
                v = OScalars -> GetTuple1(GetId(x,Dim[1]-dy+y,z,Dim));
                BScalars -> SetTuple1(GetId(x+dx,y,z+dz,NewDim),v);
            }
        }
    }
    for (x=0;x<Dim[0];x++) {
        for (y=0;y<Dim[1];y++) {
            for (z=0;z<dz;z++) {
                v = OScalars -> GetTuple1(GetId(x,y,z,Dim));
                BScalars -> SetTuple1(GetId(x+dx,y+dy,z+dz+Dim[2],NewDim),v);
                v = OScalars -> GetTuple1(GetId(x,y,Dim[2]-dz+z,Dim));
                BScalars -> SetTuple1(GetId(x+dx,y+dy,z,NewDim),v);

            }
        }
    }
   
    for (id=0;id<N;id++) {
        v = OScalars -> GetTuple1(id);
        BScalars -> SetTuple1(GetId(dx+GetX(id,Dim),dy+GetY(id,Dim),dz+GetZ(id,Dim),NewDim),v);
    }

    BScalars -> Modified();

    vtkImageData *ImageWBorder = vtkImageData::New();
    ImageWBorder -> GetPointData() -> SetScalars(BScalars);
    ImageWBorder -> SetDimensions(NewDim);

    #ifdef DEBUG
        printf("Border added!\n");
    #endif

    return ImageWBorder;
}

vtkImageData *RemoveBorder(int *DimB, vtkDoubleArray *BScalars, int *DimO) {
    #ifdef DEBUG
        printf("Removinf Border...\n");
    #endif

    double v;
    unsigned long int id;
    int dx = int(0.1*DimO[0]);
    int dy = int(0.1*DimO[1]);
    int dz = int(0.1*DimO[2]);
    unsigned long int N = DimO[0]*DimO[1]*DimO[2];
    vtkDoubleArray *OScalars = vtkDoubleArray::New();
    OScalars -> SetNumberOfTuples(N);
    for (id=0;id<N;id++) {
        v = BScalars -> GetTuple1(GetId(dx+GetX(id,DimO),dy+GetY(id,DimO),dz+GetZ(id,DimO),DimB));
        OScalars -> SetTuple1(id,v);
    }
    OScalars -> Modified();
    
    vtkImageData *Image = vtkImageData::New();
    Image -> GetPointData() -> SetScalars(OScalars);
    Image -> SetDimensions(DimO);
    #if (VTK_MAJOR_VERSION==5)
        Image -> Update();
    #endif

    #ifdef DEBUG
        printf("Border Removed!\n");
    #endif

    return Image;
}

/* ================================================================
   ROUTINES FOR FOURIER TRANSFORM APPROCH
=================================================================*/

void GetGaussianDerivativeKernel(int derivative, vtkImageData *Kernel, double sigma) {
    #ifdef DEBUG
        printf("Calculating Gaussian Derivatives (Fourier)...\n");
    #endif

    double x, y, z, d2, v;
    double s2 = pow(sigma,2);
    unsigned long int idr;
    unsigned long int register id;
    int *Dim = Kernel -> GetDimensions();
    unsigned long int N = Kernel -> GetNumberOfPoints();
    vtkDataArray *Values = Kernel -> GetPointData() -> GetScalars();
    for ( id = 0; id < N; id++ ) {
        x = GetX(id,Dim); y = GetY(id,Dim); z = GetZ(id,Dim);
        idr = GetReflectedId(x,y,z,Dim);
        x -= 0.5*Dim[0]; y -= 0.5*Dim[1]; z -= 0.5*Dim[2];
        d2 = pow(x,2) + pow(y,2) + pow(z,2);
        switch (derivative) {
            case 1: v = (x*x-s2) / pow(sigma,2.5); break;
            case 2: v = (y*y-s2) / pow(sigma,2.5); break;
            case 3: v = (z*z-s2) / pow(sigma,2.5); break;
            case 4: v = x*y / pow(sigma,4.0); break;
            case 5: v = x*z / pow(sigma,4.0); break;
            case 6: v = y*z / pow(sigma,4.0); break;
            default:v = 1.0; break;
        }
        v *= exp(-0.5*d2/s2) / pow(2*3.1415*s2,1.5);
        Values -> SetTuple1(idr,v*pow(sigma,0.25)); // Lindeberg Scaling Factor
    }
    Values -> Modified();
}

void GetImageDerivativeFourier(int derivative, double sigma, vtkImageData *ImageData, vtkDoubleArray *D) {
    #ifdef DEBUG
        printf("Calculating Image Derivatives (Fourier)...\n");
    #endif

    vtkSmartPointer<vtkImageCast> Cast = vtkSmartPointer<vtkImageCast>::New();
    #if (VTK_MAJOR_VERSION==5)
        Cast -> SetInput(ImageData);
    #else
        Cast -> SetInputData(ImageData);
    #endif
    Cast -> SetOutputScalarTypeToDouble();
    Cast -> Update();

    vtkImageData *Gauss = Cast -> GetOutput();
    GetGaussianDerivativeKernel(derivative,Gauss,sigma); 

    #ifdef DEBUG
        printf("Calculating Fourier Transform...\n");
    #endif

    vtkSmartPointer<vtkImageFFT> FFTImage = vtkSmartPointer<vtkImageFFT>::New();
    FFTImage -> SetDimensionality(3);
    #if (VTK_MAJOR_VERSION==5)
        FFTImage -> SetInput(ImageData);
    #else
        FFTImage -> SetInputData(ImageData);
    #endif
    FFTImage -> Update();    
    vtkImageData *ImageF = FFTImage -> GetOutput();

    vtkSmartPointer<vtkImageFFT> FFTGauss = vtkSmartPointer<vtkImageFFT>::New();
    FFTGauss -> SetDimensionality(3);
    #if (VTK_MAJOR_VERSION==5)
        FFTGauss -> SetInput(Gauss);
    #else
        FFTGauss -> SetInputData(Gauss);
    #endif
    FFTGauss -> Update();    
    vtkImageData *GaussF = FFTGauss -> GetOutput();

    vtkSmartPointer<vtkImageData> Mult = vtkSmartPointer<vtkImageData>::New();
    Mult -> DeepCopy(ImageF);

    double re1, re2, im1, im2;
    unsigned long int register id;

    for (id=0;id<ImageF->GetNumberOfPoints();id++) {
        re1 = ImageF -> GetPointData() -> GetScalars() -> GetComponent(id,0);
        im1 = ImageF -> GetPointData() -> GetScalars() -> GetComponent(id,1);
        re2 = GaussF -> GetPointData() -> GetScalars() -> GetComponent(id,0);
        im2 = GaussF -> GetPointData() -> GetScalars() -> GetComponent(id,1);
        Mult -> GetPointData() -> GetScalars() -> SetComponent(id,0,re1*re2-im1*im2);
        Mult -> GetPointData() -> GetScalars() -> SetComponent(id,1,re1*im2+re2*im1);
    }

    vtkSmartPointer<vtkImageRFFT> RFFT = vtkSmartPointer<vtkImageRFFT>::New();
    #if (VTK_MAJOR_VERSION==5)
        RFFT -> SetInput(Mult);
    #else
        RFFT -> SetInputData(Mult);
    #endif
    RFFT -> Update();

    vtkSmartPointer<vtkImageExtractComponents> real = vtkSmartPointer<vtkImageExtractComponents>::New();
    #if (VTK_MAJOR_VERSION==5)
        real -> SetInput(RFFT -> GetOutput());
    #else
        real -> SetInputData(RFFT -> GetOutput());
    #endif
    real -> SetComponents(0);
    real -> Update();

    D -> DeepCopy(real->GetOutput()->GetPointData()->GetScalars());
    D -> Modified();
    return;
}

void GetHessianEigenvaluesFourier(double sigma, vtkImageData *Image, vtkDoubleArray *L1, vtkDoubleArray *L2, vtkDoubleArray *L3) {
    #ifdef DEBUG
        printf("Calculating Hessian Eigenvalues (Fourier)...\n");
    #endif

    unsigned long int id;
    unsigned long int N = Image -> GetNumberOfPoints();
    double H[3][3], Eva[3], Eve[3][3], dxx, dyy, dzz, dxy, dxz, dyz, l1, l2, l3, frobnorm;
    vtkSmartPointer<vtkDoubleArray> Dxx = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> Dyy = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> Dzz = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> Dxy = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> Dxz = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> Dyz = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> Fro = vtkSmartPointer<vtkDoubleArray>::New();
    Dxx -> SetNumberOfTuples(N);
    Dyy -> SetNumberOfTuples(N);
    Dzz -> SetNumberOfTuples(N);
    Dxy -> SetNumberOfTuples(N);
    Dxz -> SetNumberOfTuples(N);
    Dyz -> SetNumberOfTuples(N);
    Fro -> SetNumberOfTuples(N);
    GetImageDerivativeFourier(1,sigma,Image,Dxx);
    GetImageDerivativeFourier(2,sigma,Image,Dyy);
    GetImageDerivativeFourier(3,sigma,Image,Dzz);
    GetImageDerivativeFourier(4,sigma,Image,Dxy);
    GetImageDerivativeFourier(5,sigma,Image,Dxz);
    GetImageDerivativeFourier(6,sigma,Image,Dyz);
    for ( id = N; id--; ) {
        l1 = l2 = l3 = 0.0;
        H[0][0]=Dxx->GetTuple1(id); H[0][1]=Dxy->GetTuple1(id); H[0][2]=Dxz->GetTuple1(id);
        H[1][0]=Dxy->GetTuple1(id); H[1][1]=Dyy->GetTuple1(id); H[1][2]=Dyz->GetTuple1(id);
        H[2][0]=Dxz->GetTuple1(id); H[2][1]=Dyz->GetTuple1(id); H[2][2]=Dzz->GetTuple1(id);
        frobnorm = FrobeniusNorm(H);
        if (H[0][0]+H[1][1]+H[2][2]<0.0) {
            vtkMath::Diagonalize3x3(H,Eva,Eve);
            l1 = Eva[0]; l2 = Eva[1]; l3 = Eva[2];
            Sort(&l1,&l2,&l3);
        }
        L1 -> SetTuple1(id,l1);
        L2 -> SetTuple1(id,l2);
        L3 -> SetTuple1(id,l3);
        Fro -> SetTuple1(id,frobnorm);
    }
    double ftresh,frobenius_norm_range[2];
    Fro -> GetRange(frobenius_norm_range);
    ftresh = sqrt(frobenius_norm_range[1]);
    for ( id = N; id--; ) {
        if ( Fro->GetTuple1(id) < ftresh) {
            L1 -> SetTuple1(id,0.0);
            L2 -> SetTuple1(id,0.0);
            L3 -> SetTuple1(id,0.0);
        }
    }
    L1 -> Modified();
    L2 -> Modified();
    L3 -> Modified();
}

/* ================================================================
   ROUTINES FOR DISCRETE APPROCH
=================================================================*/

void GetImageDerivativeDiscrete(vtkDataArray *Image, int *dim, char direction, vtkFloatArray *Derivative) {
    #ifdef DEBUG
        printf("Calculating Image Derivatives (Discrete)...\n");
    #endif

    int i, j, k;
    double d, f1, f2;
    if (direction=='x') {
        for (i = dim[0]; i--;) {
            for (j = dim[1]; j--;) {
                for (k = dim[2]; k--;) {
                    if (i==0) {
                        Image -> GetTuple(GetId(1,j,k,dim),&f1);
                        Image -> GetTuple(GetId(0,j,k,dim),&f2);
                    }
                    if (i==dim[0]-1) {
                        Image -> GetTuple(GetId(dim[0]-1,j,k,dim),&f1);
                        Image -> GetTuple(GetId(dim[0]-2,j,k,dim),&f2);
                    }
                    if (i>0&i<dim[0]-1) {
                        Image -> GetTuple(GetId(i+1,j,k,dim),&f1);
                        Image -> GetTuple(GetId(i-1,j,k,dim),&f2);
                        f1 /= 2.0; f2 /= 2.0;
                    }
                    d = f1 - f2;
                    Derivative -> SetTuple(GetId(i,j,k,dim),&d);
                }
            }
        }
    }
    if (direction=='y') {
        for (i = dim[0]; i--;) {
            for (j = dim[1]; j--;) {
                for (k = dim[2]; k--;) {
                    if (j==0) {
                        Image -> GetTuple(GetId(i,1,k,dim),&f1);
                        Image -> GetTuple(GetId(i,0,k,dim),&f2);
                    }
                    if (j==dim[1]-1) {
                        Image -> GetTuple(GetId(i,dim[1]-1,k,dim),&f1);
                        Image -> GetTuple(GetId(i,dim[1]-2,k,dim),&f2);
                    }
                    if (j>0&j<dim[1]-1) {
                        Image -> GetTuple(GetId(i,j+1,k,dim),&f1);
                        Image -> GetTuple(GetId(i,j-1,k,dim),&f2);
                        f1 /= 2.0; f2 /= 2.0;
                    }
                    d = f1 - f2;
                    Derivative -> SetTuple(GetId(i,j,k,dim),&d);
                }
            }
        }
    }
    if (direction=='z') {
        for (i = dim[0]; i--;) {
            for (j = dim[1]; j--;) {
                for (k = dim[2]; k--;) {
                    if (k==0) {
                        Image -> GetTuple(GetId(i,j,1,dim),&f1);
                        Image -> GetTuple(GetId(i,j,0,dim),&f2);
                    }
                    if (k==dim[2]-1) {
                        Image -> GetTuple(GetId(i,j,dim[2]-1,dim),&f1);
                        Image -> GetTuple(GetId(i,j,dim[2]-2,dim),&f2);
                    }
                    if (k>0&k<dim[2]-1) {
                        Image -> GetTuple(GetId(i,j,k+1,dim),&f1);
                        Image -> GetTuple(GetId(i,j,k-1,dim),&f2);
                        f1 /= 2.0; f2 /= 2.0;
                    }
                    d = f1 - f2;
                    Derivative -> SetTuple(GetId(i,j,k,dim),&d);
                }
            }
        }
    }
}

void GetHessianEigenvaluesDiscrete(double sigma, vtkImageData *Image, vtkDoubleArray *L1, vtkDoubleArray *L2, vtkDoubleArray *L3) {
    #ifdef DEBUG
        printf("Calculating Hessian Eigeinvalues (Discrete)...\n");
    #endif

    int *Dim = Image -> GetDimensions();
    unsigned long int id;
    unsigned long int N = Image -> GetNumberOfPoints();
    double H[3][3], Eva[3], Eve[3][3], dxx, dyy, dzz, dxy, dxz, dyz, l1, l2, l3, frobnorm;

    #ifdef DEBUG
        printf("Calculating Gaussian Convolution...\n");
    #endif

    vtkSmartPointer<vtkImageGaussianSmooth> Gauss = vtkSmartPointer<vtkImageGaussianSmooth>::New();
    #if (VTK_MAJOR_VERSION==5)
        Gauss -> SetInput(Image);
    #else
        Gauss -> SetInputData(Image);
    #endif
    Gauss -> SetDimensionality(3);
    Gauss -> SetRadiusFactors(10,10,10);
    Gauss -> SetStandardDeviations(sigma,sigma,sigma);
    Gauss -> Update();

    vtkFloatArray *ImageG = (vtkFloatArray*) Gauss -> GetOutput() -> GetPointData() -> GetScalars();

    vtkSmartPointer<vtkFloatArray> Dx = vtkSmartPointer<vtkFloatArray>::New(); Dx -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Dy = vtkSmartPointer<vtkFloatArray>::New(); Dy -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Dz = vtkSmartPointer<vtkFloatArray>::New(); Dz -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Dxx = vtkSmartPointer<vtkFloatArray>::New(); Dxx -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Dyy = vtkSmartPointer<vtkFloatArray>::New(); Dyy -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Dzz = vtkSmartPointer<vtkFloatArray>::New(); Dzz -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Dxy = vtkSmartPointer<vtkFloatArray>::New(); Dxy -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Dxz = vtkSmartPointer<vtkFloatArray>::New(); Dxz -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Dyz = vtkSmartPointer<vtkFloatArray>::New(); Dyz -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Fro = vtkSmartPointer<vtkFloatArray>::New();
    Fro -> SetNumberOfComponents(1);
    Fro -> SetNumberOfTuples(N);

    GetImageDerivativeDiscrete(ImageG,Dim,'x',Dx);
    GetImageDerivativeDiscrete(ImageG,Dim,'y',Dy);
    GetImageDerivativeDiscrete(ImageG,Dim,'z',Dz);
    GetImageDerivativeDiscrete(Dx,Dim,'x',Dxx);
    GetImageDerivativeDiscrete(Dy,Dim,'y',Dyy);
    GetImageDerivativeDiscrete(Dz,Dim,'z',Dzz);
    GetImageDerivativeDiscrete(Dy,Dim,'x',Dxy);
    GetImageDerivativeDiscrete(Dz,Dim,'x',Dxz);
    GetImageDerivativeDiscrete(Dz,Dim,'y',Dyz);

    for ( id = N; id--; ) {
        l1 = l2 = l3 = 0.0;
        H[0][0]=Dxx->GetTuple1(id); H[0][1]=Dxy->GetTuple1(id); H[0][2]=Dxz->GetTuple1(id);
        H[1][0]=Dxy->GetTuple1(id); H[1][1]=Dyy->GetTuple1(id); H[1][2]=Dyz->GetTuple1(id);
        H[2][0]=Dxz->GetTuple1(id); H[2][1]=Dyz->GetTuple1(id); H[2][2]=Dzz->GetTuple1(id);
        frobnorm = FrobeniusNorm(H);
        if (H[0][0]+H[1][1]+H[2][2]<0.0) {
            vtkMath::Diagonalize3x3(H,Eva,Eve);
            l1 = Eva[0]; l2 = Eva[1]; l3 = Eva[2];
            Sort(&l1,&l2,&l3);
        }
        L1 -> SetTuple1(id,l1);
        L2 -> SetTuple1(id,l2);
        L3 -> SetTuple1(id,l3);
        Fro -> SetTuple1(id,frobnorm);
    }
    double ftresh,frobenius_norm_range[2];
    Fro -> GetRange(frobenius_norm_range);
    ftresh = sqrt(frobenius_norm_range[1]);
    for ( id = N; id--; ) {
        if ( Fro->GetTuple1(id) < ftresh) {
            L1 -> SetTuple1(id,0.0);
            L2 -> SetTuple1(id,0.0);
            L3 -> SetTuple1(id,0.0);
        }
    }
    L1 -> Modified();
    L2 -> Modified();
    L3 -> Modified();
}

/* ================================================================
   VESSELNESS ROUTINE
=================================================================*/

void GetVesselness(double sigma, vtkImageData *Image, vtkDoubleArray *L1, vtkDoubleArray *L2, vtkDoubleArray *L3) {

    #ifdef DEBUG
        printf("Calculating Vesselness...\n");
    #endif

    double c = 500.0;
    double beta = 0.5;
    double alpha = 0.5;
    double std = 2 * c * c;
    double rbd = 2 * beta * beta;
    double rad = 2 * alpha * alpha;
    double l1, l2, l3, ra, ran, rb, rbn, st, stn, ft_old, ft_new;
    unsigned long int id;
    unsigned long int N = Image -> GetNumberOfPoints();

    //GetHessianEigenvaluesFourier(sigma,Image,L1,L2,L3);
    GetHessianEigenvaluesDiscrete(sigma,Image,L1,L2,L3);
    
    for ( id = N; id--; ) {
        l1 = L1 -> GetTuple1(id);
        l2 = L2 -> GetTuple1(id);
        l3 = L3 -> GetTuple1(id);
        if (l2<0&&l3<0) {

            ra = fabs(l2) / fabs(l3);
            ran = -ra * ra;

            rb = fabs(l1) / sqrt(l2*l3);
            rbn = -rb * rb;

            st = sqrt(l1*l1+l2*l2+l3*l3);
            stn = -st * st;

            ft_new = (1-exp(ran/rad)) * exp(rbn/rbd) * (1-exp(stn/std));

            //L1 is used to return vesselness values
            L1 -> SetTuple1(id,ft_new);

        } else L1 -> SetTuple1(id,0.0);
    }
    L1 -> Modified();
}

/* ================================================================
   DIVERGENT FILTER
=================================================================*/

void GetDivergentFilter(int *Dim, vtkDoubleArray *Scalars) {

    #ifdef DEBUG
        printf("Calculating Divergent Filter...\n");
    #endif

    int register j, i;
    int x, y, z, s = 2;
    unsigned long int id;
    double v, norm, V[6][3];
    int Dx[6] = {1,-1,0,0,0,0};
    int Dy[6] = {0,0,1,-1,0,0};
    int Dz[6] = {0,0,0,0,1,-1};
    int MI[3][3] = {{1,0,0},{0,1,0},{0,0,1}};

    vtkSmartPointer<vtkDoubleArray> Div = vtkSmartPointer<vtkDoubleArray>::New();
    Div -> SetNumberOfTuples(Scalars->GetNumberOfTuples());

    for (z = s+1; z < Dim[2]-s-1; z++) {
        for (y = s+1; y < Dim[1]-s-1; y++) {
            for (x = s+1; x < Dim[0]-s-1; x++) {
                v = 0;
                id = GetId(x,y,z,Dim);
                if (Scalars->GetTuple1(id)) {
                    for (i = 0; i < 6; i++) {
                        for (j = 0; j < 3; j++) {
                            V[i][j]  = Scalars -> GetTuple1(GetId(x+s*Dx[i]+MI[j][0],y+s*Dy[i]+MI[j][1],z+s*Dz[i]+MI[j][2],Dim));
                            V[i][j] -= Scalars -> GetTuple1(GetId(x+s*Dx[i]-MI[j][0],y+s*Dy[i]-MI[j][1],z+s*Dz[i]-MI[j][2],Dim));
                        }
                        norm = sqrt(pow(V[i][0],2)+pow(V[i][1],2)+pow(V[i][2],2));
                        if (norm) {V[i][0]/=norm; V[i][1]/=norm; V[i][2]/=norm; }
                    }
                    v = (V[0][0]-V[1][0])+(V[2][1]-V[3][1])+(V[4][2]-V[5][2]);
                    v = (v<0) ? -v / 6.0 : 0.0;
                }
                Div -> InsertTuple1(id,v);
            }
        }
    }
    Div -> Modified();
    Scalars -> DeepCopy(Div);
    Scalars -> Modified();

}

/* ================================================================
   MAIN ROUTINE
=================================================================*/

int main(int argc, char *argv[]) {     
           
    printf("===============\n");
    printf("---MitoGraph---\n");
    printf("===============\n");
    printf("File Name: %s\n",argv[1]);
    printf("Scales to run: [%1.3f-%1.3f]\n",1.00,1.50);
    printf("==============\n\n");

    // Loading multi-paged TIFF file (Supported by VTK 6.2 and higher)
    char _prefix[64];
    sprintf(_prefix,"%s.tif",argv[1]);
    vtkSmartPointer<vtkTIFFReader> TIFFReader = vtkSmartPointer<vtkTIFFReader>::New();
    TIFFReader -> SetFileName(_prefix);
    TIFFReader -> Update();

    //vtkSmartPointer<vtkStructuredPointsReader> reader = vtkSmartPointer<vtkStructuredPointsReader>::New();
    //reader -> SetFileName(argv[1]);
    //reader -> Update();
    //vtkImageData *ImageO = reader -> GetOutput();

    vtkImageData *ImageO = Convert16To8bit(TIFFReader->GetOutput());

    if (!ImageO) printf("Format not supported.\n");

    vtkImageData *Image = AddBorder(ImageO);

    int N = Image -> GetNumberOfPoints();

    vtkSmartPointer<vtkDoubleArray> AUX1 = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> AUX2 = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> AUX3 = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> VSSS = vtkSmartPointer<vtkDoubleArray>::New();

    AUX1 -> SetNumberOfTuples(N);
    AUX2 -> SetNumberOfTuples(N);
    AUX3 -> SetNumberOfTuples(N);
    VSSS -> SetNumberOfTuples(N);
    VSSS -> FillComponent(0,0);

    unsigned long int id;
    double sigma, vn, vo;
    int *Dim = Image -> GetDimensions();

    for ( sigma = 1.00; sigma <= 1.60; sigma += 0.10 ) {
        
        printf("Sigma = %1.3f\n",sigma);
        
        GetVesselness(sigma,Image,AUX1,AUX2,AUX3);
        
        for ( id = N; id--; ) {
            vn = AUX1 -> GetTuple1(id);
            vo = VSSS -> GetTuple1(id);
            if ( vn > vo ) {
                VSSS -> SetTuple1(id,vn);
            }
        }

    }
    VSSS -> Modified();

    GetDivergentFilter(Dim,VSSS);

    int *DimO = ImageO -> GetDimensions();

    vtkImageData *ImageEnhanced = RemoveBorder(Dim,VSSS,DimO);

    SaveImageData(ImageEnhanced,"ImageData_Vesselness.vtk");

    vtkSmartPointer<vtkContourFilter> Filter = vtkSmartPointer<vtkContourFilter>::New();
    #if (VTK_MAJOR_VERSION==5)
        Filter -> SetInput(ImageEnhanced);
    #else
        Filter -> SetInputData(ImageEnhanced);
    #endif
    Filter -> SetValue(1,0.16667);
    Filter -> Update();

    SavePolyData(Filter->GetOutput(),"Surface.vtk");

    return 0;
}

