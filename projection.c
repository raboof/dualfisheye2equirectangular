/*
MIT License

Copyright (c) 2016 F.P. Sluiter

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <getopt.h>
#include <libgen.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
So the issue is:
- source picture and target picture should have different sizes
- source picture and target picture can have different centeres, so a
translation is in order.
- Lets start with:
- target picture is in phi/theta, It should have size pi,2pi. But cuts should be
possible (lets consider that, phi theta does not allow for pretty cuts, does
it).
*/

/* Compile with: gcc projection.c -Wall -o project -lm
 * math.h does not like to be linked directly...
 */

/* Flag set by ‘--verbose’. */
static int verbose_flag;

typedef struct double2 {
  double x;
  double y;
} double2;
typedef struct double3 {
  double x;
  double y;
  double z;
} double3;

enum CameraMode { FRONT, FRONT235, BACK, UP, DOWN, SAMSUNG_GEAR_360, THETAS };

typedef struct configuration {
  char *xmap_filename;
  char *ymap_filename;
  int xmap_set;
  int ymap_set;
  int rows;   // target
  int cols;   // target
  int height; // source
  int width;  // source
  int rows_set;
  int cols_set;
  int height_set;
  int width_set;
  int crop;
  enum CameraMode mode;
} configuration;

configuration parse_options(int argc, char **argv) {
  int c;
  configuration po; // to hold parsed options
  po.xmap_filename = NULL;
  po.ymap_filename = NULL;
  po.xmap_set = 0;
  po.ymap_set = 0;
  po.rows = 0;
  po.cols = 0;
  po.rows_set = 0;
  po.cols_set = 0;
  po.height_set = 0;
  po.width_set = 0;
  po.mode = FRONT; // default

  while (1) {
    static struct option long_options[] = {
        /* These options set a flag. */
        {"verbose", no_argument, &verbose_flag, 1},
        {"brief", no_argument, &verbose_flag, 0},

        /* These options don’t set a flag.
           We distinguish them by their indices. */
        {"help", no_argument, 0, 'q'},

        /* options with arg*/
        {"xmap", required_argument, 0, 'x'},
        {"ymap", required_argument, 0, 'y'},
        {"rows", required_argument, 0, 'r'},   // target
        {"cols", required_argument, 0, 'c'},   // target
        {"height", required_argument, 0, 'h'}, // source
        {"width", required_argument, 0, 'w'},  // source
        {"mode", required_argument, 0, 'm'},
        {"crop", required_argument, 0, 'b'},

        {0, 0, 0, 0}};
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long(argc, argv, "qx:y:r:c:h:w:m:b:", long_options,
                    &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c) {
    case 0:
      /* If this option set a flag, do nothing else now. */
      if (long_options[option_index].flag != 0)
        break;
      printf("option %s", long_options[option_index].name);
      if (optarg)
        printf(" with arg %s", optarg);
      printf("\n");
      break;

    case 'x':
      po.xmap_filename = optarg;
      po.xmap_set++;
      break;

    case 'y':
      po.ymap_filename = optarg;
      po.ymap_set++;
      break;
    case 'h':
      po.height = atoi(optarg);
      po.height_set++;
      break;
    case 'w':
      po.width = atoi(optarg);
      po.width_set++;
      break;
    case 'c':
      po.cols = atoi(optarg);
      po.cols_set++;
      break;
    case 'r':
      po.rows = atoi(optarg);
      po.rows_set++;
      break;
    case 'b':
      po.crop = atoi(optarg);
      break;
    case 'm':
      if (strcmp(optarg, "front") == 0) {
        po.mode = FRONT;
      } else if (strcmp(optarg, "front235") == 0) {
        po.mode = FRONT235;
      } else if (strcmp(optarg, "down") == 0) {
        po.mode = DOWN;
      } else if (strcmp(optarg, "samsung_gear_360") == 0) {
        po.mode = SAMSUNG_GEAR_360;
      } else if (strcmp(optarg, "thetas") == 0) {
        po.mode = THETAS;
      }
      /* more else if clauses */
      else /* default: */
      {
        printf("Camera mode %s not implemented \n", optarg);
        exit(1);
      }
      break;

    case '?':
    /* getopt_long already printed an error message. */
    case 'q':
      printf("Usage: %s -x|--xmap FILE_x.pgm -y|--ymap FILE_y.pgm -h|--height "
             "300 -w|--width 400 -r|--rows 600 -c|--cols 800 \n",
             argv[0]);
      printf("h,w is source size, r,c is targetsize \n");
      exit(1);
      break;

    default:
      abort();
    }
  }

  /* Instead of reporting ‘--verbose’
     and ‘--brief’ as they are encountered,
     we report the final status resulting from them. */
  if (verbose_flag) {
    switch (po.mode) {
    case FRONT:
      printf("Camera: Front proj\n");
      break;
    case FRONT235:
      printf("Camera: Front 235 proj\n");
      break;
    case DOWN:
      printf("Camera: Down proj\n");
      break;
    case SAMSUNG_GEAR_360:
      printf("Camera: samsung_gear_360\n");
      break;
    case THETAS:
      printf("Camera: Theta S\n");
      break;
    default:
      printf("Camera mode not in verbose, exiting\n");
      exit(1);
    }
  }
  /* Print any remaining command line arguments (not options). */
  if (optind < argc) {
    printf("ERROR: non-option ARGV-elements: ");
    while (optind < argc)
      printf("%s ", argv[optind++]);
    putchar('\n');
    exit(1);
  }

  if (po.xmap_set != 1 || po.ymap_set != 1) {
    printf("ERROR: Xmap and ymap are mandatory arguments and have to appear "
           "only once!\ntry --help for help\n\n ");
    exit(-1);
  }
  if (po.rows_set != 1 || po.cols_set != 1) {
    printf("ERROR: Target Rows and Cols are mandatory arguments and have to "
           "appear only once!\ntry --help for help\n\n ");
    exit(-1);
  }
  if (po.height_set != 1 || po.width_set != 1) {
    printf("ERROR: Source Height and Width are mandatory arguments and have to "
           "appear only once!\ntry --help for help\n\n ");
    exit(-1);
  }
  return po;
}
#define MAXROWS 9000
#define MAXCOLS 9000

int pgmWrite_ASCII(char *filename, int rows, int cols, int **image,
                   char *comment_string) {
  FILE *file; /* pointer to the file buffer */
              //      int maxval;        /* maximum value in the image array */
  long nwritten = 0; /* counter for the number of pixels written */
  long x, y;         /* for loop counters */

  /* return 0 if the dimensions are larger than the image array. */
  if (rows > MAXROWS || cols > MAXCOLS) {
    printf("ERROR: row/col specifications larger than image array:\n");
    return (0);
  }

  /* open the file; write header and comments specified by the user. */
  if ((file = fopen(filename, "w")) == NULL) {
    printf("ERROR: file open failed\n");
    return (0);
  }
  fprintf(file, "P2\n");

  if (comment_string != NULL)
    fprintf(file, "# %s \n", comment_string);

  /* write the dimensions of the image */
  fprintf(file, "%i %i \n", cols, rows);

  /* NOTE: MAXIMUM VALUE IS WHITE; COLOURS ARE SCALED FROM 0 - */
  /* MAXVALUE IN A .PGM FILE. */

  /* WRITE MAXIMUM VALUE TO FILE */
  fprintf(file, "%d\n", (int)65535);

  /* Write data */

  for (y = 0; y < rows; y++) {
    for (x = 0; x < cols; x++) {
      fprintf(file, "%i ", image[y][x]);
      nwritten++;
    }
    fprintf(file, "\n");
  }
  fprintf(file, "\n");

  printf("\nNumber of pixels total(from rows*cols): %i\n", rows * cols);
  printf("Number of pixels written in file %s: %ld\n\n", filename, nwritten);

  fclose(file);
  return (1);
}

/* So, to get the x’,y’ position for the circular image we will have to first
 * pass the
 * coordinates x,y from the rectangular output image to spherical coordinates
 * using the
 * first coordinate system, then those to the second shown spherical coordinate
 * system,
 * then those to the polar projection and then pass the polar system to cardinal
 * x’,y’.
 */
double2 evaluatePixel_Front(double2 outCoord, double2 srcSize) {
  double2 o;
  double theta, phi;
  double3 sphericCoords;
  double phi2_over_pi;
  double theta2;
  double2 inCentered;

  // convert outcoords to radians (180 = pi, so half a sphere)
  o.x = outCoord.x / srcSize.x;
  o.y = outCoord.y / srcSize.y;
  theta = (1.0 - o.x) * M_PI;
  phi = o.y * M_PI;

  // Convert outcoords to spherical (x,y,z on unisphere)
  sphericCoords.x = cos(theta) * sin(phi);
  sphericCoords.y = sin(theta) * sin(phi);
  sphericCoords.z = cos(phi);

  // Convert spherical to input coordinates...
  theta2 = atan2(-sphericCoords.z, sphericCoords.x);
  phi2_over_pi = acos(sphericCoords.y) / (M_PI);

  inCentered.x = ((phi2_over_pi * cos(theta2)) + 0.5) * srcSize.x;
  inCentered.y = ((phi2_over_pi * sin(theta2)) + 0.5) * srcSize.y;

  return inCentered;
}

void gen_front_maps(configuration cfg, int **image_x, int **image_y) {
  int x, y;

  printf("Front proj\n");
  for (y = 0; y < cfg.rows; y++) {
    for (x = 0; x < cfg.cols; x++) {
      double2 o = evaluatePixel_Front(
          (double2){
              ((double)x / ((double)cfg.cols)) * ((cfg.width) - (2 * cfg.crop)),
              ((double)y / (double)cfg.rows) * ((cfg.height) - (2 * cfg.crop))},
          (double2){(cfg.width) - (2 * cfg.crop), cfg.height - (2 * cfg.crop)});
      image_x[y][x] = (int)round(o.x) + cfg.crop;
      image_y[y][x] = (int)round(o.y) + cfg.crop;
    }
  }
}

void gen_samsung_gear_360_maps(configuration cfg, int **image_x,
                               int **image_y) {
  int x, y;
  // two circels, but with overlap/crop
  // idear: rescale the x,y square. So circel should be smaller (2*crop) than
  // actual image size.
  // so srcSize in 2 dimensions is: - 2*crop
  // and translate with 1* crop for left and 2* crop for right
  printf("samsung_gear_360 with crop %i\n", cfg.crop);

  for (y = 0; y < cfg.rows; y++) {
    for (x = 0; x < cfg.cols / 2; x++) {

      double2 o =
          evaluatePixel_Front((double2){((double)x / ((double)cfg.cols / 2)) *
                                            ((cfg.width / 2) - (2 * cfg.crop)),
                                        ((double)y / (double)cfg.rows) *
                                            ((cfg.height) - (2 * cfg.crop))},
                              (double2){(cfg.width / 2) - (2 * cfg.crop),
                                        cfg.height - (2 * cfg.crop)});
      if (o.x < 0)
        o.x = 0;
      if (o.y < 0)
        o.y = 0;
      if (o.x >= cfg.width / 2)
        o.x = (cfg.width / 2) - 1;
      if (o.y >= cfg.height)
        o.y = cfg.height - 1;

      // crop on outsides
      // for samsung gear 360 it should be ~38
      image_x[y][x] = (int)round(o.x) + cfg.crop;
      image_x[y][x + cfg.cols / 2] = (int)round(o.x) + cfg.crop + cfg.width / 2;
      image_y[y][x] = (int)round(o.y) + cfg.crop;
      image_y[y][x + cfg.cols / 2] = (int)round(o.y) + cfg.crop;
    }
  }
}

void gen_thetas_maps(configuration cfg, int **image_x, int **image_y) {
  int x, y;
  // two circels, but with overlap/crop
  // idear: rescale the x,y square. So circel should be smaller (2*crop) than
  // actual image size.
  // so srcSize in 2 dimensions is: - 2*crop
  // and translate with 1* crop for left and 2* crop for right
  printf("ThetaS with crop %i\n", cfg.crop);

  if (cfg.rows != cfg.cols / 2) {
    printf("ERROR: ThetaS projection only works if cols == 2*rows \ntry --help "
           "for help\n\n ");
    exit(-1);
  }

  for (y = 0; y < cfg.rows; y++) {
    for (x = 0; x < cfg.cols / 2; x++) {

      double2 o =
          evaluatePixel_Front((double2){((double)x / ((double)cfg.cols / 2)) *
                                            ((cfg.width / 2) - (2 * cfg.crop)),
                                        ((double)y / (double)cfg.rows) *
                                            ((cfg.height) - (2 * cfg.crop))},
                              (double2){(cfg.width / 2) - (2 * cfg.crop),
                                        cfg.height - (2 * cfg.crop)});
      if (o.x < 0)
        o.x = 0;
      if (o.y < 0)
        o.y = 0;
      if (o.x >= cfg.width / 2)
        o.x = (cfg.width / 2) - 1;
      if (o.y >= cfg.height)
        o.y = cfg.height - 1;

      // crop on outsides
      // for thetas it should be ~55
      /*
      left side:
      Y goes from 0 to 959 and becomes X 959 to 0: (for y=0;y<=959;y++) x =
      959-y;
      X goes from 0 to 959 and becomes Y 0 to 959: (for x = 0; x<=959;x++) y =x;
      right side:
      Y goes from 0 to 959 and becomes X 960 to 1919 (for y =0; y<=959;y++) x =
      960+y;
      X goes from 960 to 1919 and becomes Y 959 to 0 (for x = 960; x <= 1919;
      x++) y = 1919-x
      */

      // Left:
      image_x[y][(cfg.height - 1) - x] = (int)round(o.y) + cfg.crop;
      image_y[y][(cfg.height - 1) - x] = (int)round(o.x) + cfg.crop;
      // Right:
      image_x[(cfg.height - 1) - y][x + cfg.cols / 2] =
          (int)round(o.y) + cfg.crop + cfg.width / 2;
      image_y[(cfg.height - 1) - y][x + cfg.cols / 2] =
          (int)round(o.x) + cfg.crop;
    }
  }
}

/*
double3 rotate(double3 p, double3 o, double angle){
  double3 r;

  r.x = [u(ux + vy + wz)(1-cos th) + x cos th + (-wy + vz)sin th]
  r.y = [v(ux + vy + wz)(1-cos th) + y cos th + ( wy - uz)sin th]
  r.z = [w(ux + vy + wz)(1-cos th) + z cos th + (-vy + uy)sin th]


}
*/
/* So, to get the x’,y’ position for the circular image we will have to first
 * pass the
 * coordinates x,y from the rectangular output image to spherical coordinates
 * using the
 * first coordinate system, then those to the second shown spherical coordinate
 * system,
 * then those to the polar projection and then pass the polar system to cardinal
 * x’,y’.
 */
double3 evaluatePixel_down(double2 outCoord, double2 srcSize, double2 trgSize) {
  double2 o;
  double theta, phi;
  double3 sphericCoords;
  double phi2_over_pi;
  double theta2;
  double3 inCentered;

  o.x = outCoord.x / trgSize.x;
  o.y = outCoord.y / (trgSize.y - 1);

  theta = (1 - o.x) * 2 * M_PI; // quadrant??
  phi = o.y * M_PI;

  sphericCoords.x = cos(theta) * sin(phi);
  sphericCoords.y = sin(theta) * sin(phi);
  sphericCoords.z = cos(phi);
  // rotated over x sphericoords are here
  // x2 = x
  // y2 = -z
  // z2 = y
  // that was done on below old, where they use atan2 in stead of?
  // old
  /*	theta2 = atan2(-sphericCoords.z,sphericCoords.x);
          phi2_over_pi = acos(sphericCoords.y)/(M_PI);
  */
  theta2 = atan2(-sphericCoords.z, sphericCoords.x);
  phi2_over_pi = acos(sphericCoords.y) /
                 (M_PI); // klopt dit nog? delen door pi?? gaat om y?

  inCentered.x =
      ((phi2_over_pi * cos(theta2)) + 0.5) * srcSize.x; //+0.5 voor middelpunt?
  inCentered.y =
      ((phi2_over_pi * sin(theta2)) + 0.5) * srcSize.y; //+0.5 voor middelpunt?
  inCentered.z = phi2_over_pi;

  // dst = sample(src,inCentered);
  printf(" [%i %i: (%f %f) (%f %f %f) (%f %f %f) ] \n", (int)floor(outCoord.x),
         (int)floor(outCoord.y), o.x, o.y, sphericCoords.x, sphericCoords.y,
         sphericCoords.z, inCentered.x, inCentered.y, inCentered.z);
  return inCentered;
}

void gen_down_maps(configuration cfg, int **image_x, int **image_y) {
  int x, y;

  printf("Front proj\n");
  for (y = 0; y < cfg.rows; y++) {
    for (x = 0; x < cfg.cols; x++) {
      double3 o =
          evaluatePixel_down((double2){((double)x), ((double)y)},
                             (double2){(cfg.width), cfg.height},
                             (double2){(double)cfg.cols, (double)cfg.rows});
      //      if(y>cfg.rows/2){
      image_x[y][x] = (int)round(o.x);
      image_y[y][x] = (int)round(o.y);
      //      } else {
      //          image_x[y][x] = 0;
      //          image_y[y][x] = 0;
      //      }
    }
  }
}

/* This fron 235 function to accomodate kodaks wider fisheye.
 * So, to get the x’,y’ position for the circular image we will have to first
 * pass the
 * coordinates x,y from the rectangular output image to spherical coordinates
 * using the
 * first coordinate system, then those to the second shown spherical coordinate
 * system,
 * then those to the polar projection and then pass the polar system to cardinal
 * x’,y’.
 */
double2 evaluatePixel_Front235(double2 outCoord, double2 srcSize) {
  double2 o;
  double theta, phi;
  double3 sphericCoords;
  double phi2_over_pi;
  double theta2;
  double2 inCentered;

  // convert outcoords to radians (180 = pi, so half a sphere)
  o.x = outCoord.x / srcSize.x;
  o.y = outCoord.y / srcSize.y;
  theta = (1.0 - o.x) * M_PI * 1.3055; // 235/180=1.3055
  phi = o.y * M_PI * 1.3055;           // 235/180=1.3055

  // Convert outcoords to spherical (x,y,z on unisphere)
  sphericCoords.x = cos(theta) * sin(phi);
  sphericCoords.y = sin(theta) * sin(phi);
  sphericCoords.z = cos(phi);

  // Convert spherical to input coordinates...
  theta2 = atan2(-sphericCoords.z, sphericCoords.x);
  phi2_over_pi = acos(sphericCoords.y) / (M_PI * 1.3055);

  inCentered.x = ((phi2_over_pi * cos(theta2)) + 0.5) * srcSize.x;
  inCentered.y = ((phi2_over_pi * sin(theta2)) + 0.5) * srcSize.y;

  return inCentered;
}

void gen_front235_maps(configuration cfg, int **image_x, int **image_y) {
  int x, y;

  printf("Front 235 proj\n");
  for (y = 0; y < cfg.rows; y++) {
    for (x = 0; x < cfg.cols; x++) {
      double2 o = evaluatePixel_Front(
          (double2){
              ((double)x / ((double)cfg.cols)) * ((cfg.width) - (2 * cfg.crop)),
              ((double)y / (double)cfg.rows) * ((cfg.height) - (2 * cfg.crop))},
          (double2){(cfg.width) - (2 * cfg.crop), cfg.height - (2 * cfg.crop)});
      image_x[y][x] = (int)round(o.x) + cfg.crop;
      image_y[y][x] = (int)round(o.y) + cfg.crop;
    }
  }
}

int main(int argc, char **argv) {
  int y;
  int **image_x;
  int **image_y;
  configuration cfg = parse_options(argc, argv);

  if (cfg.xmap_filename)
    printf("xmapfile: %s\n", cfg.xmap_filename);
  if (cfg.ymap_filename)
    printf("ymapfile: %s\n", cfg.ymap_filename);

  image_x = malloc((cfg.rows) * sizeof(*image_x));
  for (y = 0; y < (cfg.rows); y++)
    image_x[y] = malloc((cfg.cols) * sizeof(*(image_x[y])));
  image_y = malloc((cfg.rows) * sizeof(*image_y));
  for (y = 0; y < (cfg.rows); y++)
    image_y[y] = malloc((cfg.cols) * sizeof(*(image_y[y])));

  switch (cfg.mode) {
  case FRONT:
    gen_front_maps(cfg, image_x, image_y);
    break;
  case FRONT235:
    gen_front235_maps(cfg, image_x, image_y);
    break;
  case DOWN:
    gen_down_maps(cfg, image_x, image_y);
    break;
  case SAMSUNG_GEAR_360:
    gen_samsung_gear_360_maps(cfg, image_x, image_y);
    break;
  case THETAS:
    gen_thetas_maps(cfg, image_x, image_y);
    break;
  default:
    printf("Camera mode not implemented\n");
    exit(1);
  }

  printf("Writing files\n");
  pgmWrite_ASCII(cfg.ymap_filename, cfg.rows, cfg.cols, image_y,
                 cfg.ymap_filename);
  pgmWrite_ASCII(cfg.xmap_filename, cfg.rows, cfg.cols, image_x,
                 cfg.xmap_filename);

  if (image_y) {
    for (y = 0; y < cfg.rows; y++)
      free(image_y[y]);
    free(image_y);
  }
  if (image_x) {
    for (y = 0; y < cfg.rows; y++)
      free(image_x[y]);
    free(image_x);
  }

  exit(0);
}
