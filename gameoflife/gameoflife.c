#include <endian.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <setjmp.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <sys/time.h>


#include <omp.h>

#define calcIndex(width, x,y)  ((y)*(width) + (x))

long TimeSteps = 100;

void writeParallelVTK(long timestep, int w, int h, int px, int py)
{
	int nx = w / px; 
	int ny = h / py;
	
  char filename[2048];

  snprintf(filename, sizeof(filename), "%s-%05ld%s", "parallel", timestep, ".pvti");
  FILE *fp = fopen(filename, "w");

  fprintf(fp, "<?xml version=\"1.0\"?>\n");
  fprintf(fp, "<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  fprintf(fp, "<PImageData WholeExtent=\"%d %d %d %d %d %d\" Origin=\"0 0 0\" Spacing=\"%le %le %le\">\n", 0, w, 0, h, 0, 0, 1.0, 1.0, 0.0);
  fprintf(fp, "<PCellData Scalars=\"%s\">\n", "gol");
  fprintf(fp, "<PDataArray type=\"Float32\" Name=\"%s\" format=\"appended\" offset=\"0\"/>\n", "gol");
  fprintf(fp, "</PCellData>\n");
  for (int i = 0; i < nx * ny; i++)
  {
    int xOffset = (i % nx) * px;
	int yOffset = (i / nx) * py;

    fprintf(fp, "<Piece Extent=\"%d %d %d %d %d %d\" Source=\"%s-%d-%05ld%s\"/>", xOffset, xOffset + px, yOffset, yOffset + py, 0, 0, "gol", i, timestep, ".vti");
  }

  fprintf(fp, "</PImageData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}

void writeVTK2(long timestep, double *data, char prefix[1024], int w, int h, int Gw, int offsetX, int offsetY, int num) {
  char filename[2048];  
  int x,y; 
  
  //int offsetX=0;
  //int offsetY=0;
  float deltax=1.0;
  long  nxy = w * h * sizeof(float);  

  snprintf(filename, sizeof(filename), "%s-%d-%05ld%s", prefix, num, timestep, ".vti");
  FILE* fp = fopen(filename, "w");

  fprintf(fp, "<?xml version=\"1.0\"?>\n");
  fprintf(fp, "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  fprintf(fp, "<ImageData WholeExtent=\"%d %d %d %d %d %d\" Origin=\"0 0 0\" Spacing=\"%le %le %le\">\n", offsetX, offsetX + w, offsetY, offsetY + h, 0, 0, deltax, deltax, 0.0);
  fprintf(fp, "<CellData Scalars=\"%s\">\n", prefix);
  fprintf(fp, "<DataArray type=\"Float32\" Name=\"%s\" format=\"appended\" offset=\"0\"/>\n", prefix);
  fprintf(fp, "</CellData>\n");
  fprintf(fp, "</ImageData>\n");
  fprintf(fp, "<AppendedData encoding=\"raw\">\n");
  fprintf(fp, "_");
  fwrite((unsigned char*)&nxy, sizeof(long), 1, fp);

  for (y = offsetY; y < h + offsetY; y++) {
    for (x = offsetX; x < w + offsetX; x++) {
      float value = data[calcIndex(Gw, x,y)];
      fwrite((unsigned char*)&value, sizeof(float), 1, fp);
    }
  }
  
  fprintf(fp, "\n</AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}


void show(double* currentfield, int w, int h) {
  printf("\033[H");
  int x,y;
  for (y = 0; y < h; y++) {
    for (x = 0; x < w; x++) printf(currentfield[calcIndex(w, x,y)] ? "\033[07m  \033[m" : "  ");
    //printf("\033[E");
    printf("\n");
  }
  fflush(stdout);
}

void evolve(double* currentfield, double* newfield, int w, int h, int pX, int pY, long t) {
  writeParallelVTK(t, w, h, pX, pY);
  #pragma omp parallel for collapse(2) //schedule(static, 1) 
  for(int rectangleX = 0; rectangleX < w / pX; rectangleX++) {
    for(int rectangleY = 0; rectangleY < h / pY; rectangleY++) {
      int offsetX = pX * rectangleX;
      int offsetY = pY * rectangleY;               
      for (int y = offsetY; y < offsetY + pY; y++) {
        for (int x = offsetX; x <  offsetX + pX; x++) {
          int neighbourCount = 0;

          for (int y1 = y - 1; y1 <= y + 1; y1++)
            for (int x1 = x - 1; x1 <= x + 1; x1++)
              if (x1 >= 0 && x1 < w && y1 >= 0 && y1 < h && (x1 != x || y1 != y) && currentfield[calcIndex(w, x1, y1)])
                neighbourCount++;

          newfield[calcIndex(w, x,y)] = neighbourCount == 3 || currentfield[calcIndex(w, x,y)] && neighbourCount == 2;
        }
      }
      writeVTK2(t,currentfield,"gol", pX, pY, w, offsetX, offsetY, (rectangleY * (w / pX)) + rectangleX);
    }
  }
}
 
void filling(double* currentfield, int w, int h, char *path) {
  if (path == 0){
    int i;
    for (i = 0; i < h*w; i++) {
      currentfield[i] = (rand() < RAND_MAX / 10) ? 1 : 0; ///< init domain randomly
    }
    return;
  }
  FILE *file = fopen(path, "r");
  if(file == 0){
    printf("File not found");
    exit(0);
  }

  while(getc(file) == '#'){
    while(getc(file) != '\n');
  }

  char c;
  ///x
  while(getc(file) != '=');
  getc(file);
  int x = 0;
  while((c = getc(file)) != ','){
    x = x * 10 + c - '0';
  }

  ///y
  while(getc(file) != '=');
  getc(file);
  int y = 0;
  while((c = getc(file)) != ','){
    y = y * 10 + c - '0';
  }

  if(x != w || y != h){
    fclose(file);
    printf("Wrong sizes");
    exit(0);
  }

  while(getc(file) != '\n');

  x = 0;
  y = 0;
  c = getc(file);
  while(c != '!'){
    while(c != '$' && c != '!'){
      int rl = 1;
      if(c >= '0' && c <= '9'){
        rl = 0;
        while(c != 'b' && c != 'o'){
          rl = rl * 10 + c - '0';
          c = getc(file);    
          while(c == '\n'){
            c = getc(file);
          }
        }
      } 
      for(int index = 0; index < rl; index++) {
        currentfield[calcIndex(w, x, y)] = c == 'o';
        x++;
      }
      c = getc(file);
      if(c == '\n'){
        c = getc(file);
      }
    }
    if(c == '$'){
      c = getc(file);
      if(c == '\n'){
        c = getc(file);
      }
    }
    x = 0;
    y++;
  }
}
 
void game(int w, int h, int pX, int pY, char *path) {
  double *currentfield = calloc(w*h, sizeof(double));
  double *newfield     = calloc(w*h, sizeof(double));
  
  //printf("size unsigned %d, size long %d\n",sizeof(float), sizeof(long));
  
  filling(currentfield, w, h, path);
  long t;
  for (t=0;t<TimeSteps;t++) {
    
    //show(currentfield, w, h); 
    
    evolve(currentfield, newfield, w, h, pX, pY, t);
    
    //printf("%ld timestep\n",t);
    //writeVTK2(t,currentfield,"gol", w, h);
    
    //usleep(200000);

    //SWAP
    double *temp = currentfield;
    currentfield = newfield;
    newfield = temp;
  }
  
  free(currentfield);
  free(newfield);
  
}
 
int main(int c, char **v) {
  int nX = 0, nY = 0, pX = 0, pY = 0;
  char *path = 0;
  if (c > 1) TimeSteps = atoi(v[1]); ///< read timesteps
  if (c > 2) nX = atoi(v[2]); ///< read nX
  if (c > 3) nY = atoi(v[3]); ///< read nY

  if (c > 4) pX = atoi(v[4]); ///< read pX
  if (c > 5) pY = atoi(v[5]); ///< read pY

  if (c > 6) path = v[6]; ///< read path
  
  if (nX <= 0) nX = 2; ///< default nX
  if (nY <= 0) nY = 2; ///< default nY
  if (pX <= 0) pX = 8; ///< default width
  if (pY <= 0) pY = 8; ///< default height

  int w = nX * pX;
  int h = nY * pY;

  //printf("Start Game of Life with w=%d h=%d\n", w, h);

  game(w, h, pX, pY, path);
}