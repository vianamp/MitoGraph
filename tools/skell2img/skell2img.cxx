/*

    Export skeleton edges with length larger than threshold as a binary image

*/

#include "includes.h"

int main(int argc, char *argv[]) {     

    int i;
    std::string _skpath;
    std::string _rfpath;
    double _dxy, _dz, _length_threshold;

    // Collecting input parameters
    for (i = 0; i < argc; i++) {
        if (!strcmp(argv[i],"-skell")) {
            _skpath = argv[i+1];
        }
        if (!strcmp(argv[i],"-reference")) {
            _rfpath = argv[i+1];
        }
        if (!strcmp(argv[i],"-threshold")) {
            _length_threshold = (double)atof(argv[i+1]);
        }
        if (!strcmp(argv[i],"-xy")) {
            _dxy = atof(argv[i+1]);
        }
        if (!strcmp(argv[i],"-z")) {
            _dz = atof(argv[i+1]);
        }
    }

    printf("Skeleton: %s\n",_skpath.c_str());
    printf("Reference Image: %s\n",_rfpath.c_str());
    printf("Threshold: %1.3f\n",_length_threshold);

    //
    // Load Reference Image
    //

    vtkSmartPointer<vtkTIFFReader> TIFFReader = vtkSmartPointer<vtkTIFFReader>::New();
    TIFFReader -> SetFileName((_rfpath+".tif").c_str());
    TIFFReader -> Update();
    
    vtkSmartPointer<vtkImageData> Image = TIFFReader -> GetOutput();

    Image -> GetPointData() -> GetScalars() -> FillComponent(0,0);

    //
    // Load Skeleton
    //

    vtkSmartPointer<vtkPolyDataReader> PolyReader = vtkSmartPointer<vtkPolyDataReader>::New();
    PolyReader -> SetFileName(_skpath.c_str());
    PolyReader -> Update();

    vtkSmartPointer<vtkPolyData> Skell = PolyReader -> GetOutput();
    
    printf("Number of Polydata Points: %d\n",(int)Skell->GetNumberOfPoints());

    vtkIdType id, id2;
    double r[3], length;
    vtkSmartPointer<vtkCell> Edge;

    for (id = 0; id < Skell->GetNumberOfPoints(); id++) {
        Skell -> GetPoint(id,r);
        length = Skell -> GetPointData() -> GetArray("Length") -> GetTuple1(id);

        if (length > _length_threshold) {
            Image -> GetPointData() -> GetScalars() -> SetTuple1(Image->FindPoint((int)r[0],(int)r[1],0), 255);
        }

    }

    vtkSmartPointer<vtkTIFFWriter> Writer = vtkSmartPointer<vtkTIFFWriter>::New();
    Writer -> SetInputData(Image);
    Writer -> SetFileName((_rfpath+"-mask.tif").c_str());
    Writer -> Write();

    return 0;
}
