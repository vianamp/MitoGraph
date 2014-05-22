#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vtkStructuredPointsReader.h>
#include <vtkImageData.h>
#include <vtkImageFFT.h>
#include <vtkDataObject.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredPoints.h>
#include <vtkDataArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkContourFilter.h>
#include <vtkPointData.h>
#include <vtkPolyDataWriter.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkStructuredPointsWriter.h>
//#include <vtkImageFourierCenter.h>
#include <vtkImageRFFT.h>
#include <vtkImageCast.h>
#include <vtkMath.h>
#include <vtkImageExtractComponents.h>

//#define DEBUG

int  GetX(unsigned long int id, int *Dim);
int  GetY(unsigned long int id, int *Dim);
int  GetZ(unsigned long int id, int *Dim);
unsigned long int GetId(int x, int y, int z, int *Dim);
void fill_gaussian_derivative(int derivative, vtkImageData *Kernel, double sigma);
void get_image_derivative(int derivative, vtkImageData *Image, vtkDoubleArray *D);
int get_reflected_id(int x, int y, int z, int *Dim);
void get_Hessian_eigenvalues(double sigma, vtkImageData *Image, vtkDoubleArray *L1, vtkDoubleArray *L2, vtkDoubleArray *L3);
void swap(double *x, double *y);
void sort(double *l1, double *l2, double *l3);
void get_vesselness(double sigma, vtkImageData *Image, vtkDoubleArray *L1, vtkDoubleArray *L2, vtkDoubleArray *L3);

/**========================================================
 Auxiliar functions
 =========================================================*/

int  GetX(unsigned long int id, int *Dim)              { return (int) id%Dim[0]; }
int  GetY(unsigned long int id, int *Dim)              { return (int)(id%(Dim[0]*Dim[1]))/Dim[0]; }
int  GetZ(unsigned long int id, int *Dim)              { return (int) id/(Dim[0]*Dim[1]); }
unsigned long int GetId(int x, int y, int z, int *Dim) { return       x+y*Dim[0]+z*Dim[0]*Dim[1]; }

void IMAGEDATA_DIVERGENT_FILTER(int *Dim, vtkDoubleArray *Scalars) {

    int register j, i;
    int x, y, z, id, s = 2;
    double v, norm, V[6][3];
    int Dx[6] = {1,-1,0,0,0,0};
    int Dy[6] = {0,0,1,-1,0,0};
    int Dz[6] = {0,0,0,0,1,-1};
    int MI[3][3] = {{1,0,0},{0,1,0},{0,0,1}};

    vtkSmartPointer<vtkDoubleArray> Div = vtkSmartPointer<vtkDoubleArray>::New();
    Div -> SetNumberOfTuples(Scalars->GetNumberOfTuples());

    double max=0;
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
                max = (v>max) ? v : max;
                Div -> InsertTuple1(id,v);
            }
        }
    }
    Div -> Modified();
    Scalars -> DeepCopy(Div);
    Scalars -> Modified();
    printf("MAX = %f\n",max);
}

void _save_imagedata(vtkImageData *Image) {
#ifdef DEBUG
    printf("Saving ImageData File...\n");
#endif
    vtkSmartPointer<vtkStructuredPointsWriter> writer = vtkSmartPointer<vtkStructuredPointsWriter>::New();
    writer->SetInputData(Image);
    writer->SetFileName("temp00.vtk");
    writer->Write();
#ifdef DEBUG
    printf("File Saved!\n");
#endif
}

vtkImageData *remove_border(int *DimB, vtkDoubleArray *BScalars, int *DimO) {
    int id;
    double v;
    int dx = int(0.1*DimO[0]);
    int dy = int(0.1*DimO[1]);
    int dz = int(0.1*DimO[2]);
    int N = DimO[0]*DimO[1]*DimO[2];
    vtkDoubleArray *OScalars = vtkDoubleArray::New();
    OScalars -> SetNumberOfTuples(N);
    for (id=0;id<N;id++) {
        v = BScalars -> GetTuple1(GetId(dx+GetX(id,DimO),dy+GetY(id,DimO),dz+GetZ(id,DimO),DimB));
        OScalars -> SetTuple1(id,v);
    }
    OScalars -> Modified();
    
    vtkImageData *Temp = vtkImageData::New();
    //Temp -> SetScalarTypeToDouble();
    Temp -> GetPointData() -> SetScalars(OScalars);
    Temp -> SetDimensions(DimO);
    //Temp -> Update();

    return Temp;
}

vtkImageData *add_border(vtkImageData *ImageData) {
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

    printf("%d\t%d\t%d\n",Lx,Ly,Lz);

    vtkDataArray *OScalars = ImageData -> GetPointData() -> GetScalars();
    vtkDoubleArray *BScalars = vtkDoubleArray::New();
    BScalars -> SetNumberOfTuples(Lx*Ly*Lz);
    BScalars -> FillComponent(0,0);

    double v;
    int x, y, z, id;
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

#ifdef DEBUG
    printf("Creating New ImageData...\n");
#endif

    vtkImageData *ImageWBorder = vtkImageData::New();
    //ImageWBorder -> SetScalarTypeToDouble();
    ImageWBorder -> GetPointData() -> SetScalars(BScalars);
    ImageWBorder -> SetDimensions(NewDim);

    //_save_imagedata(ImageWBorder);

    return ImageWBorder;
}

int get_reflected_id(int x, int y, int z, int *Dim) {
    int rx = ceil(0.5*Dim[0]);
    int ry = ceil(0.5*Dim[1]);
    int rz = ceil(0.5*Dim[2]);
    int sx = (x-(rx-0.5)<0) ? -rx : Dim[0]-rx;
    int sy = (y-(ry-0.5)<0) ? -ry : Dim[1]-ry;
    int sz = (z-(rz-0.5)<0) ? -rz : Dim[2]-rz;
    return GetId(x-sx,y-sy,z-sz,Dim);
}

void fill_xygaussian(vtkImageData *ImageData, double sigma) {
    int ir;
    int register i;
    double x, y, z, d2, v;
    double s2 = pow(sigma,2);
    int *Dim = ImageData -> GetDimensions();
    int N = ImageData -> GetNumberOfPoints();
    vtkDataArray *Values = ImageData -> GetPointData() -> GetScalars();
    for ( i = 0; i < N; i++ ) {
        x = GetX(i,Dim); y = GetY(i,Dim); z = GetZ(i,Dim);
        x -= 0.5*Dim[0]; y -= 0.5*Dim[1]; z -= 0.5*Dim[2];
        d2 = pow(x,2) + pow(y,2) + pow(z,2);
        v = 255 * exp(-0.5*(x*x+y*y)/s2);
        Values -> SetTuple1(i,(int)v); // Lindeberg Scaling Factor
    }
    Values -> Modified();
}

void fill_xzgaussian(vtkImageData *ImageData, double sigma) {
    int ir;
    int register i;
    double x, y, z, d2, v;
    double s2 = pow(sigma,2);
    int *Dim = ImageData -> GetDimensions();
    int N = ImageData -> GetNumberOfPoints();
    vtkDataArray *Values = ImageData -> GetPointData() -> GetScalars();
    for ( i = 0; i < N; i++ ) {
        x = GetX(i,Dim); y = GetY(i,Dim); z = GetZ(i,Dim);
        x -= 0.5*Dim[0]; y -= 0.5*Dim[1]; z -= 0.5*Dim[2];
        d2 = pow(x,2) + pow(y,2) + pow(z,2);
        v = 255 * exp(-0.5*(x*x+z*z)/s2);
        Values -> SetTuple1(i,(int)v); // Lindeberg Scaling Factor
    }
    Values -> Modified();
}


void fill_gaussian_derivative(int derivative, vtkImageData *Kernel, double sigma) {
    int ir;
    int register i;
    double x, y, z, d2, v;
    double s2 = pow(sigma,2);
    int *Dim = Kernel -> GetDimensions();
    int N = Kernel -> GetNumberOfPoints();
    vtkDataArray *Values = Kernel -> GetPointData() -> GetScalars();
    for ( i = 0; i < N; i++ ) {
        x = GetX(i,Dim); y = GetY(i,Dim); z = GetZ(i,Dim);
        ir = get_reflected_id(x,y,z,Dim);
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
        Values -> SetTuple1(ir,v*pow(sigma,0.25)); // Lindeberg Scaling Factor
    }
    Values -> Modified();
}

void get_image_derivative(int derivative, double sigma, vtkImageData *ImageData, vtkDoubleArray *D) {
    vtkSmartPointer<vtkImageCast> Cast = vtkSmartPointer<vtkImageCast>::New();
    Cast -> SetInputData(ImageData);
    Cast -> SetOutputScalarTypeToDouble();
    Cast -> Update();

    vtkImageData *Gauss = Cast -> GetOutput();
    fill_gaussian_derivative(derivative,Gauss,sigma); 

    vtkSmartPointer<vtkImageFFT> FFTImage = vtkSmartPointer<vtkImageFFT>::New();
    FFTImage -> SetDimensionality(3);
    FFTImage -> SetInputData(ImageData);
    FFTImage -> Update();    
    vtkImageData *ImageF = FFTImage -> GetOutput();

    vtkSmartPointer<vtkImageFFT> FFTGauss = vtkSmartPointer<vtkImageFFT>::New();
    FFTGauss -> SetDimensionality(3);
    FFTGauss -> SetInputData(Gauss);
    FFTGauss -> Update();    
    vtkImageData *GaussF = FFTGauss -> GetOutput();

    vtkSmartPointer<vtkImageData> Mult = vtkSmartPointer<vtkImageData>::New();
    Mult -> DeepCopy(ImageF);

    int register i;
    double re1, re2, im1, im2;
    // Try to use the complex multiplication from VTK
    for (i=0;i<ImageF->GetNumberOfPoints();i++) {
        re1 = ImageF -> GetPointData() -> GetScalars() -> GetComponent(i,0);
        im1 = ImageF -> GetPointData() -> GetScalars() -> GetComponent(i,1);
        re2 = GaussF -> GetPointData() -> GetScalars() -> GetComponent(i,0);
        im2 = GaussF -> GetPointData() -> GetScalars() -> GetComponent(i,1);
        Mult -> GetPointData() -> GetScalars() -> SetComponent(i,0,re1*re2-im1*im2);
        Mult -> GetPointData() -> GetScalars() -> SetComponent(i,1,re1*im2+re2*im1);
    }

    vtkSmartPointer<vtkImageRFFT> RFFT = vtkSmartPointer<vtkImageRFFT>::New();
    RFFT -> SetInputData(Mult);
    RFFT -> Update();

    vtkSmartPointer<vtkImageExtractComponents> real = vtkSmartPointer<vtkImageExtractComponents>::New();
    real -> SetInputData(RFFT -> GetOutput());
    real -> SetComponents(0);
    real -> Update();

    D -> DeepCopy(real->GetOutput()->GetPointData()->GetScalars());
    D -> Modified();
    return;
}

double frobenius_norm(double M[3][3]) {
    double f = 0.0;
    for (int i = 3;i--;)
        for (int j = 3;j--;)
            f += M[i][j]*M[i][j];
    return sqrt(f);
}

void get_Hessian_eigenvalues_Fourier(double sigma, vtkImageData *Image, vtkDoubleArray *L1, vtkDoubleArray *L2, vtkDoubleArray *L3) {
    int i, N = Image -> GetNumberOfPoints();
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
    get_image_derivative(1,sigma,Image,Dxx);
    get_image_derivative(2,sigma,Image,Dyy);
    get_image_derivative(3,sigma,Image,Dzz);
    get_image_derivative(4,sigma,Image,Dxy);
    get_image_derivative(5,sigma,Image,Dxz);
    get_image_derivative(6,sigma,Image,Dyz);
    for ( i = N; i--; ) {
        l1 = l2 = l3 = 0.0;
        H[0][0]=Dxx->GetTuple1(i); H[0][1]=Dxy->GetTuple1(i); H[0][2]=Dxz->GetTuple1(i);
        H[1][0]=Dxy->GetTuple1(i); H[1][1]=Dyy->GetTuple1(i); H[1][2]=Dyz->GetTuple1(i);
        H[2][0]=Dxz->GetTuple1(i); H[2][1]=Dyz->GetTuple1(i); H[2][2]=Dzz->GetTuple1(i);
        frobnorm = frobenius_norm(H);
        if (H[0][0]+H[1][1]+H[2][2]<0.0) {
            vtkMath::Diagonalize3x3(H,Eva,Eve);
            l1 = Eva[0]; l2 = Eva[1]; l3 = Eva[2];
            sort(&l1,&l2,&l3);
        }
        L1 -> SetTuple1(i,l1);
        L2 -> SetTuple1(i,l2);
        L3 -> SetTuple1(i,l3);
        Fro -> SetTuple1(i,frobnorm);
    }
    double ftresh,frobenius_norm_range[2];
    Fro -> GetRange(frobenius_norm_range);
    ftresh = sqrt(frobenius_norm_range[1]);
    for ( i = N; i--; ) {
        if ( Fro->GetTuple1(i) < ftresh) {
            L1 -> SetTuple1(i,0.0);
            L2 -> SetTuple1(i,0.0);
            L3 -> SetTuple1(i,0.0);
        }
    }
    L1 -> Modified();
    L2 -> Modified();
    L3 -> Modified();
}

//==========================================================================

void get_image_derivative_discrete(vtkDataArray *Image, int *dim, char direction, vtkFloatArray *Derivative) {
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

void get_Hessian_eigenvalues_discrete(double sigma, vtkImageData *Image, vtkDoubleArray *L1, vtkDoubleArray *L2, vtkDoubleArray *L3) {

    int *Dim = Image -> GetDimensions();
    int i, N = Image -> GetNumberOfPoints();
    double H[3][3], Eva[3], Eve[3][3], dxx, dyy, dzz, dxy, dxz, dyz, l1, l2, l3, frobnorm;

#ifdef DEBUG
    printf("Gaussian Convolution...\n");
#endif

    vtkSmartPointer<vtkImageGaussianSmooth> Gauss = vtkSmartPointer<vtkImageGaussianSmooth>::New();
    Gauss -> SetInputData(Image);
    Gauss -> SetDimensionality(3);
    Gauss -> SetRadiusFactors(10,10,10);
    Gauss -> SetStandardDeviations(sigma,sigma,sigma);
    Gauss -> Update();

#ifdef DEBUG
    printf("Gaussian Convolution Done!\n");
#endif

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

    get_image_derivative_discrete(ImageG,Dim,'x',Dx);
    get_image_derivative_discrete(ImageG,Dim,'y',Dy);
    get_image_derivative_discrete(ImageG,Dim,'z',Dz);
    get_image_derivative_discrete(Dx,Dim,'x',Dxx);
    get_image_derivative_discrete(Dy,Dim,'y',Dyy);
    get_image_derivative_discrete(Dz,Dim,'z',Dzz);
    get_image_derivative_discrete(Dy,Dim,'x',Dxy);
    get_image_derivative_discrete(Dz,Dim,'x',Dxz);
    get_image_derivative_discrete(Dz,Dim,'y',Dyz);

    for ( i = N; i--; ) {
        l1 = l2 = l3 = 0.0;
        H[0][0]=Dxx->GetTuple1(i); H[0][1]=Dxy->GetTuple1(i); H[0][2]=Dxz->GetTuple1(i);
        H[1][0]=Dxy->GetTuple1(i); H[1][1]=Dyy->GetTuple1(i); H[1][2]=Dyz->GetTuple1(i);
        H[2][0]=Dxz->GetTuple1(i); H[2][1]=Dyz->GetTuple1(i); H[2][2]=Dzz->GetTuple1(i);
        frobnorm = frobenius_norm(H);
        if (H[0][0]+H[1][1]+H[2][2]<0.0) {
            vtkMath::Diagonalize3x3(H,Eva,Eve);
            l1 = Eva[0]; l2 = Eva[1]; l3 = Eva[2];
            sort(&l1,&l2,&l3);
        }
        L1 -> SetTuple1(i,l1);
        L2 -> SetTuple1(i,l2);
        L3 -> SetTuple1(i,l3);
        Fro -> SetTuple1(i,frobnorm);
    }
    double ftresh,frobenius_norm_range[2];
    Fro -> GetRange(frobenius_norm_range);
    ftresh = sqrt(frobenius_norm_range[1]);
    for ( i = N; i--; ) {
        if ( Fro->GetTuple1(i) < ftresh) {
            L1 -> SetTuple1(i,0.0);
            L2 -> SetTuple1(i,0.0);
            L3 -> SetTuple1(i,0.0);
        }
    }
    L1 -> Modified();
    L2 -> Modified();
    L3 -> Modified();

}

//==========================================================================

void swap(double *x, double *y) {
    double t = *y; *y = *x; *x = t;
}
void sort(double *l1, double *l2, double *l3) {
    if (fabs(*l1) > fabs(*l2)) swap(l1,l2);
    if (fabs(*l2) > fabs(*l3)) swap(l2,l3);
    if (fabs(*l1) > fabs(*l2)) swap(l1,l2);
}

void get_vesselness(double sigma, vtkImageData *Image, vtkDoubleArray *L1, vtkDoubleArray *L2, vtkDoubleArray *L3) {

    double c = 500.0;
    double beta = 0.5;
    double alpha = 0.5;
    double std = 2 * c * c;
    double rbd = 2 * beta * beta;
    double rad = 2 * alpha * alpha;
    double l1, l2, l3, ra, ran, rb, rbn, st, stn, ft_old, ft_new;
    int i, N = Image -> GetNumberOfPoints();

#ifdef DEBUG
    printf("Calculating Hessian...\n");
#endif

    //get_Hessian_eigenvalues_Fourier(sigma,Image,L1,L2,L3);
    get_Hessian_eigenvalues_discrete(sigma,Image,L1,L2,L3);
    
#ifdef DEBUG
    printf("Hessian Done!\n");
#endif


    for ( i = N; i--; ) {
        l1 = L1 -> GetTuple1(i);
        l2 = L2 -> GetTuple1(i);
        l3 = L3 -> GetTuple1(i);
        if (l2<0&&l3<0) {

            ra = fabs(l2) / fabs(l3);
            ran = -ra * ra;

            rb = fabs(l1) / sqrt(l2*l3);
            rbn = -rb * rb;

            st = sqrt(l1*l1+l2*l2+l3*l3);
            stn = -st * st;

            ft_new = (1-exp(ran/rad)) * exp(rbn/rbd) * (1-exp(stn/std));

            L1 -> SetTuple1(i,ft_new);
        } else L1 -> SetTuple1(i,0.0);
    }
    L1 -> Modified();
}

int main(int argc, char *argv[]) {     
           
#ifdef DEBUG
    printf("File Name: %s\n",argv[1]);
#endif


    vtkSmartPointer<vtkStructuredPointsReader> reader = vtkSmartPointer<vtkStructuredPointsReader>::New();
    reader -> SetFileName(argv[1]);
    reader -> Update();

    vtkImageData *ImageO = reader -> GetOutput();

#ifdef DEBUG
    printf("Adding Border...\n");
#endif

    vtkImageData *Image = add_border(ImageO);

#ifdef DEBUG
    printf("Border Added!\n");
#endif

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

    int i;
    double sigma, vn, vo;
    int *Dim = Image -> GetDimensions();
    for (sigma=1.00;sigma<=1.60;sigma+=0.10) {
        printf("Sigma = %1.3f\n",sigma);
        get_vesselness(sigma,Image,AUX1,AUX2,AUX3);
        for ( i = N; i--; ) {
            vn = AUX1 -> GetTuple1(i);
            vo = VSSS -> GetTuple1(i);
            if (vn>vo) {
                VSSS -> SetTuple1(i,vn);
            }
        }
    }
    VSSS -> Modified();

    double tt[2];
    VSSS->GetRange(tt);
    printf("%f\t%f\n",tt[0],tt[1]);

    IMAGEDATA_DIVERGENT_FILTER(Dim,VSSS);

    VSSS->GetRange(tt);
    printf("%f\t%f\n",tt[0],tt[1]);


    int *DimO = ImageO -> GetDimensions();
    vtkSmartPointer<vtkDoubleArray> Div = vtkSmartPointer<vtkDoubleArray>::New();
    Div -> SetNumberOfTuples(ImageO->GetNumberOfPoints());

    vtkImageData *Temp = remove_border(Dim,VSSS,DimO);

    double *qq = Temp->GetScalarRange();
    printf("%f\t%f\n",qq[0],qq[1]);

    vtkSmartPointer<vtkStructuredPointsWriter> writer = vtkSmartPointer<vtkStructuredPointsWriter>::New();
    writer->SetInputData(Temp);
    writer->SetFileName("ImageData_Vesselness.vtk");
    writer->Write();


    vtkSmartPointer<vtkContourFilter> Filter = vtkSmartPointer<vtkContourFilter>::New();
    Filter -> SetInputData(Temp);
    Filter -> SetValue(1,0.16667);
    Filter -> Update();

    vtkSmartPointer<vtkPolyDataWriter> PolyDataWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
    PolyDataWriter -> SetInputData(Filter->GetOutput());
    PolyDataWriter -> SetFileName("Surface.vtk");
    PolyDataWriter -> Write();

    return 0;
}

