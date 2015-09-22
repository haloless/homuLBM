
#include <stdio.h>
#include <stdlib.h>

#include <utility>
#include <algorithm>

#define OUTPUT_CSV output_csv_
#define OUTPUT_FRONT output_front_
#define OUTPUT_VTK3D output_vtk3d_

template<typename T, size_t N>
size_t COUNTOF(const T (&arr)[N]) {
    return N;
}




template<typename T>
void writeBigEndian(FILE* fp,const T* data,size_t counts)
{
//	int width=sizeof(T);

	const char* p=(const char*)data;
	char buffer[sizeof(T)*256];
	for(size_t written=0; written<counts; written+=COUNTOF(buffer)/sizeof(T)){
		size_t toWrite=std::min(COUNTOF(buffer)/sizeof(T),counts-written);

		for(int i=0;i<toWrite;i++){
			for(int j=0;j<sizeof(T);j++){
				//swap bytes
				buffer[i*sizeof(T) + j] = p[sizeof(T)*(i+1) - j - 1];
			}
		}

		fwrite(buffer,sizeof(T) * toWrite,1,fp);
	}
}

template<typename T>
void writeBigEndianSameData(FILE* fp,const T* data,size_t counts)
{
	const char* p=(const char*)data;
	char buffer[sizeof(T)*256];
	for(int i=0;i<COUNTOF(buffer)/sizeof(T);i++){
		for(int j=0;j<sizeof(T);j++){
			//swap bytes
			buffer[i*sizeof(T) + j] = p[sizeof(T)-j-1];
		}
	}
	for(size_t i=0;i<counts;i+=COUNTOF(buffer)/sizeof(T)){
		size_t toWrite=std::min(COUNTOF(buffer)/sizeof(T),counts-i);

		fwrite(buffer,sizeof(T) * toWrite,1,fp);
	}
}

inline void writeBigEndianCastFloat(FILE* fp,const double* data,size_t count)
{
	float buffer[256];
	size_t index=0;
	for(;;){
		size_t packed=0;
		for(size_t i=0;i<COUNTOF(buffer) && i+index<count;i++,packed++){
			buffer[i]=(float)data[index+i];
		}
		writeBigEndian(fp,buffer,packed);
		if(packed<COUNTOF(buffer))break;
		index+=packed;
	}
}


static void write_cell_scalar(FILE *fp, const char *name, 
    const double *data, int nx, int ny) 
{
    fprintf(fp, "SCALARS %s float 1\n", name);
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (int j=1; j<=ny; j++) {
        for (int i=1; i<=nx; i++) {
            int idx = i + j*(nx+2);
            float f = data[idx];
            writeBigEndian(fp, &f, 1);
        }
    }
    fprintf(fp, "\n");
}
static void write_cell_scalar(FILE *fp, const char *name, 
    const double *data, int nx, int ny, int nz) 
{
    fprintf(fp, "SCALARS %s float 1\n", name);
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (int k=1; k<=nz; k++) {
        for (int j=1; j<=ny; j++) {
            for (int i=1; i<=nx; i++) {
                int idx = i + j*(nx+2) + k*(ny+2)*(nx+2);
                float f = data[idx];
                writeBigEndian(fp, &f, 1);
            }
        }
    }
    fprintf(fp, "\n");
}

extern "C" 
void OUTPUT_CSV(
    const int &step, const double &time,
    const double *color,
    const double *dens, 
    const double *visc,
    const double *pres,
    const double *uu, const double *vv,
    const double *x, const double *y, // cell position
    const double *xh, const double *yh, // face position
    const int &nx, const int &ny)
{
    char fname[256];
    sprintf(fname, "output/fluid%04d.vtk", step);
    FILE* fp = fopen(fname, "wb");
    if (!fp) {
        fprintf(stderr, "Failed to open %s\n", fname);
        exit(1);
    }
    
    // header
    fprintf(fp, "# vtk DataFile Version 2.0\n");
    // title
    fprintf(fp, "Fluid Data\n");
    // type
    fprintf(fp, "BINARY\n");
    fprintf(fp, "DATASET RECTILINEAR_GRID\n");
    // XYZ
    fprintf(fp, "DIMENSIONS %d %d %d\n", nx+1, ny+1, 1);
    fprintf(fp, "X_COORDINATES %d float\n", nx+1);
    for (int i=0; i<=nx; i++) {
        float f = xh[i];
        writeBigEndian(fp, &f, 1);
    }
    fprintf(fp, "\n");
    fprintf(fp, "Y_COORDINATES %d float\n", ny+1);
    for (int j=0; j<=ny; j++) {
        float f = yh[j];
        writeBigEndian(fp, &f, 1);
    }
    fprintf(fp, "\n");
    fprintf(fp, "Z_COORDINATES %d float\n", 1);
    {
        float f = 0.0;
        writeBigEndian(fp, &f, 1);
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "FIELD FieldData 1\n");
    fprintf(fp, "TIME 1 1 float\n");
    writeBigEndianCastFloat(fp, &time, 1);
    
    //
    fprintf(fp, "CELL_DATA %d\n", nx*ny);
    
    //write_cell_scalar(fp, "color", color, nx, ny);
    write_cell_scalar(fp, "dens", dens, nx, ny);
    //write_cell_scalar(fp, "visc", visc, nx, ny);
    //write_cell_scalar(fp, "pres", pres, nx, ny);
    
    if (0) {
    //
    fprintf(fp, "POINT_DATA %d\n", (nx+1)*(ny+1));
    
    fprintf(fp, "VECTORS vel float\n");
    for (int j=0; j<=ny; j++) {
        for (int i=0; i<=nx; i++) {
            int idx = i + j*(nx+1);
            float data[] = { uu[idx], vv[idx], 0.0 };
            writeBigEndian(fp, data, 3);
        }
    }
    fprintf(fp, "\n");
    }
    
    fclose(fp);
    
    printf("Saved %s\n", fname);
    
    return;
}

extern "C" 
void OUTPUT_FRONT(const int &step, 
    const double *xf, const double *yf,
    const double *uf, const double *vf,
    const double *tx, const double *ty,
    const int &nf)
{
    char fname[256];
    sprintf(fname, "output/front%04d.csv", step);
    
    FILE *fp = fopen(fname, "w");
    if (!fp) {
        fprintf(stderr, "Failed to open %s\n", fname);
        exit(1);
    }
    
    fprintf(fp, "x,y,z\n");
    for (int l=1; l<=(nf); l++) {
        fprintf(fp, "%lf,%lf,0.0\n", xf[l], yf[l]);
    }
    
    fclose(fp);
    printf("Saved %s\n", fname);
    return;
}

extern "C" 
void OUTPUT_VTK3D(
    const int &step, const double &time,
    const double *dens, 
    const double *xh, const double *yh, const double *zh, // face position
    const int &nx, const int &ny, const int &nz)
{
    char fname[256];
    sprintf(fname, "output/fluid%04d.vtk", step);
    FILE* fp = fopen(fname, "wb");
    if (!fp) {
        fprintf(stderr, "Failed to open %s\n", fname);
        exit(1);
    }
    
    // header
    fprintf(fp, "# vtk DataFile Version 2.0\n");
    // title
    fprintf(fp, "Fluid Data\n");
    // type
    fprintf(fp, "BINARY\n");
    fprintf(fp, "DATASET RECTILINEAR_GRID\n");
    // XYZ
    fprintf(fp, "DIMENSIONS %d %d %d\n", nx+1, ny+1, nz+1);
    fprintf(fp, "X_COORDINATES %d float\n", nx+1);
    for (int i=0; i<=nx; i++) {
        float f = xh[i];
        writeBigEndian(fp, &f, 1);
    }
    fprintf(fp, "\n");
    fprintf(fp, "Y_COORDINATES %d float\n", ny+1);
    for (int j=0; j<=ny; j++) {
        float f = yh[j];
        writeBigEndian(fp, &f, 1);
    }
    fprintf(fp, "\n");
    fprintf(fp, "Z_COORDINATES %d float\n", nz+1);
    for (int k=0; k<=nz; k++) {
        float f = zh[k];
        writeBigEndian(fp, &f, 1);
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "FIELD FieldData 1\n");
    fprintf(fp, "TIME 1 1 float\n");
    writeBigEndianCastFloat(fp, &time, 1);
    
    //
    fprintf(fp, "CELL_DATA %d\n", nx*ny*nz);
    
    //write_cell_scalar(fp, "color", color, nx, ny);
    //write_cell_scalar(fp, "dens", dens, nx, ny, nz);
    //write_cell_scalar(fp, "visc", visc, nx, ny);
    //write_cell_scalar(fp, "pres", pres, nx, ny);
    
    //
    fprintf(fp, "POINT_DATA %d\n", (nx+1)*(ny+1)*(nz+1));
    
    fprintf(fp, "SCALARS rho1 float 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (int k=0; k<=nz; k++) {
        for (int j=0; j<=ny; j++) {
            for (int i=0; i<=nx; i++) {
                float data = 0;
                for (int kk=k; kk<=k+1; kk++) {
                    for (int jj=j; jj<=j+1; jj++) {
                        for (int ii=i; ii<=i+1; ii++) {
                            int idx = ii + jj*(nx+2) + kk*(ny+2)*(nx+2);
                            data += dens[idx];
                        }
                    }
                }
                data /= 8;
                writeBigEndian(fp, &data, 1);
            }
        }
    }
    fprintf(fp, "\n");
    
    fclose(fp);
    
    printf("Saved %s\n", fname);
    
    return;
}



