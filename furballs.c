/*
 *
 *  Frolicking Furballs Safari Resort
 *  Game made by Sos within 48h
 *  Made for Ludum Dare 20 ( http://www.ludumdare.com/ )
 *
 *  some small bits of code taken from LD16 entry
 *
 *  no changes to actual game code were made after the compo
 *  except from comments and cleanup
 *
 */

////////////////////////////////////////////////////
//////////////////////////////////////////////////// INCLUDES
////////////////////////////////////////////////////

// standard includes
// included allegro to get a kickstart for the tight deadline
#include <alleggl.h>  // gl helper for allegro
#include <SOIL/SOIL.h>// Simple OpenGL Image Loader
#include <allegro.h>  // allegro - a game programming library
#include <math.h>     // just math
#include <stdio.h>    // standard

////////////////////////////////////////////////////
//////////////////////////////////////////////////// DEFINES
////////////////////////////////////////////////////

// tweakable defines
#define LOFI 1            // this disables VBOs, reduces drawing distance and furball complexity
                          // HiFi version runs smoothly on a modern PC or decent laptop
                          // LoFi version runs smoothly on a 10yo PC
#define TURN_SPEED .1     // camera/player control sensitivity
#define MOVE_SPEED 1.0    // camera/player move speed
#define WORLD_SIZE 2000.0 // side of scene square
#define GRID 256          // number of floor grid cells
#define HAYSTACK 100.0    // player height above the ground modifier

// defines for calculation simplification
#define FRAND(x) (((rand() % (int)x) * 1024.0) / 1024.0)
#define DEG(n) ((n)*180.0 / M_PI)

// player states
#define WALK 1
#define IDLE 2
#define STOP 3

#define HEIGHT 10.0      // height of camera/player
#define BOUNCE_RATE 2.0f // player jump rate
#define NUMSIZE .1f      // size of onscreen numbers
#define FUR_SLICE 8      // not used DELME

// game state
#define INTRO 1
#define PLAY 2
#define OUTRO 3

#define IQ 2       // furball bounces between change of direction
#define SPEED 1.5f // furball speed

// graphics quality settings
#if LOFI == 1
#define DRAW_DIST 500.0f // rendering / fog distance
#define PARTICLES 128    // number of max onscreen particles
#define BALLZ 666        // number of furballs
#define TRAIL 1          // hair length for hairy furballs
#define LINE_WIDTH 24.0f // line width for grass
#else
#define LINE_WIDTH 12.0f
#define PARTICLES 256
#define BALLZ 666
#define TRAIL 4
#define DRAW_DIST 2000.0f
#endif

#define MAX_ENTITIES 1024 // max scene entities
#define SAMPS 8           // number of furball voice samples
#define TALK_DELAY 100    // delay between furball sample playback

////////////////////////////////////////////////////
//////////////////////////////////////////////////// TYPEDEFS
////////////////////////////////////////////////////

typedef struct RECT RECT;
struct RECT
{
  int32_t left;
  int32_t top;
  int32_t right;
  int32_t bottom;
};

// basic vertex buffer struct
typedef struct BUFFER
{
  float *vtx, *clr;              // vertex and colour buffers
  GLuint vtx_handle, clr_handle; // VBO handles for the above
  int size,                      // number of vertices
    hardware;                    // is it VBO or array
  int frames;                    // number of frames
  int current_frame;             // currently playing frame
  GLenum mode;                   // lines or points
} BUFFER;

// blood particle struct
typedef struct BLOOD
{
  float x, y, z;    // paritcle position
  float r, g, b, a; // colour
  float vx, vy, vz; // velocity
  float s;          // size
  int alive;        // update flag
} BLOOD;

// furball (a dynamic vertex buffer + AI)
typedef struct FURBALL
{
  BUFFER *buf,               // base vertex buffer for augmentation
    *mine;                   // actually drawn buffer
  float x, y, z,             // position
    scale,                   // size
    r, g, b,                 // colour
    ultimate,                // is it dynamic (the long-haired one
    a,                       // angle of movement
    bounce,                  // bounce rate
    speed,                   // speed of movement
    bounce_rate;             // base bounce rate
  float trail[TRAIL + 1][3]; // 'hair' trail (only ultimate)
  int density,               // points per voxel
    exists,                  // render flag
    dying;                   // blood render flag
  BLOOD blood[PARTICLES];    // blood particles
  float iq,                  // directiona change interval
    smart;                   // iq base value
} FURBALL;

// a buffer with position
typedef struct ENTITY
{
  float x, y, z, s, a; // position, size, yaw (angle)
  BUFFER* buf;         // vertex buffer
} ENTITY;

////////////////////////////////////////////////////
//////////////////////////////////////////////////// GLOBALS
////////////////////////////////////////////////////

int frames = 0,   // frames rendered
  state = IDLE;   // player state
volatile int tim; // threaded timing variable

float playerx, playery, playerz, // player position
  lookx, looky, lookz, lookf;    // player look direction

BUFFER *ultimate_furball, // base buffer for ultimate furball (with ahir)
  *casual_furball;        // base buffer for simple furball

// vertex buffers for world entities
BUFFER *grass, *tall_tree, *withered_bush, *cactus, *palm,
  *stick, *pine, *hut, *fence, *church, *brickhouse;

BITMAP* ground_bmp;           // ground colour
float bounce = BOUNCE_RATE;   // player bounce rate
FURBALL* ballz[BALLZ];        // furballs
int mx, my;                   // mouse delta position
BLOOD blood[PARTICLES];       // particles
ENTITY ents[MAX_ENTITIES];    // entities (map objects) array
int num_ents = 0;             // number of entities
int game_state = INTRO;       // game state
GLuint numbers, intro, outro; // textures

SAMPLE *youknow, *furtalk[SAMPS], *die[SAMPS],
  *slash, *music, *shot, *whistle; // sounds

int talk_counter = TALK_DELAY, talking = 1; // furball talk timers
int whistle_timeout = 100;                  // furball whistle timeout

////////////////////////////////////////////////////
//////////////////////////////////////////////////// HELPER FUNCTIONS
////////////////////////////////////////////////////

// checks if GL extension is supported (stolen)
#include <string.h>
int
isExtensionSupported(const char* extension)
{
  const GLubyte* extensions = NULL;
  const GLubyte* start;
  GLubyte *where, *terminator;
  /* Extension names should not have spaces. */
  where = (GLubyte*)strchr(extension, ' ');
  if (where || *extension == '\0')
    return 0;
  extensions = glGetString(GL_EXTENSIONS);
  /* It takes a bit of care to be fool-proof about parsing the
     OpenGL extensions string. Don't be fooled by sub-strings,
     etc. */
  start = extensions;
  for (;;)
    {
      where = (GLubyte*)strstr((const char*)start, extension);
      if (!where)
        break;
      terminator = where + strlen(extension);
      if (where == start || *(where - 1) == ' ')
        if (*terminator == ' ' || *terminator == '\0')
          return 1;
      start = terminator;
    }
  return 0;
}

////////////////////////////////////////////////////
//////////////////////////////////////////////////// MATH FUNCTIONS
////////////////////////////////////////////////////

// vector length
float
length3v(float* a)
{
  return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

// 2 point distance
float
dist3v(float* a, float* b)
{
  float x, y, z;
  x = a[0] - b[0];
  y = a[1] - b[1];
  z = a[2] - b[2];
  return sqrt(x * x + y * y + z * z);
}

// dot product
float
dot3v(float* a, float* b)
{
  return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

////////////////////////////////////////////////////
//////////////////////////////////////////////////// RENDERING FUNCTIONS
////////////////////////////////////////////////////

// creates VBO from vertex array
void
vboize(BUFFER* b)
{
  if (LOFI)
    return;
  if (!isExtensionSupported("GL_ARB_vertex_buffer_object"))
    return;
  printf("Vboizing buffer with %i elements\n", b->size);
  b->hardware = 1;
  glGenBuffers(1, &(b->vtx_handle));
  glBindBuffer(GL_ARRAY_BUFFER, b->vtx_handle);
  glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 3 * b->size, b->vtx, GL_STATIC_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  glGenBuffers(1, &(b->clr_handle));
  glBindBuffer(GL_ARRAY_BUFFER, b->clr_handle);
  glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 3 * b->size, b->clr, GL_STATIC_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  printf("Removing arrays...\n");
  free(b->vtx);
  free(b->clr);
  printf("Done!\n");
}

// renders a BUFFER at desired position
void
draw_buffer_ex(BUFFER* b, float x, float y, float z, float sx, float sy, float sz, float point_size, float angle)
{
  float dist, psize, scale;
  scale = sx;
  dist = sqrt((x - playerx) * (x - playerx) + (y - playery) * (y - playery) + (z - playerz) * (z - playerz));

  // calculates point size from distance
  psize = SCREEN_W * (scale) / dist;
  if (psize < .01)
    return;
  glPointSize(psize * point_size);
  glPushMatrix();
  glScalef(sx, sy, sz);
  glTranslatef(x / sx, y / sy, z / sz);
  glRotatef(angle, 0, 1, 0);

  // draws a VBO or vertex array
  if (b->hardware)
    {
      glBindBuffer(GL_ARRAY_BUFFER, b->vtx_handle);
      glVertexPointer(3, GL_FLOAT, 0, 0);
      glBindBuffer(GL_ARRAY_BUFFER, b->clr_handle);
      glColorPointer(3, GL_FLOAT, 0, 0);
      glDrawArrays(b->mode, 0, b->size);
      glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
  else
    {
      // glBindBuffer(GL_ARRAY_BUFFER, 0);
      glVertexPointer(3, GL_FLOAT, 0, b->vtx);
      glColorPointer(3, GL_FLOAT, 0, b->clr);
      glDrawArrays(b->mode, 0, b->size);
    }
  glPopMatrix();
}

// simplified version of the above
void
draw_buffer(BUFFER* b, float x, float y, float z, float scale, float angle)
{
  draw_buffer_ex(b, x, y, z, scale, scale, scale, 1.0f, angle);
}

// draws eyes for a furball
void
draw_eyes(FURBALL* f)
{
  float x1 = 1.2, y = 10, z1 = 1, x2 = 1.2, z2 = 1, s;
  // haired furballs are bigger
  x1 *= cos(f->a) * (f->ultimate ? 1 : 10);
  z1 *= sin(f->a) * (f->ultimate ? 1 : 8);

  x2 *= cos(-f->a) * (f->ultimate ? 1 : 10);
  z2 *= sin(-f->a) * (f->ultimate ? 1 : 8);
  y = (f->ultimate ? 1.5 : 10);

  // gets point size for this furball
  glGetFloatv(GL_POINT_SIZE, &s);
  if (f->ultimate)
    {
      // gets line width if it's a hairy one
      glGetFloatv(GL_LINE_WIDTH, &s);
      s *= .25;
    }
  // if the size is too small, no use drawing eyes
  if (s < 2.0f)
    return;
  glPointSize(s * 4);
  glPushMatrix();
  glScalef(f->scale, f->scale, f->scale);
  glTranslatef(f->x / f->scale, f->y / f->scale, f->z / f->scale);

  // draws eyes
  glBegin(GL_POINTS);
  glColor3f(1, 1, 1);
  glVertex3f(x1, y, z1);
  glVertex3f(x2, y, z2);
  glEnd();
  glPointSize(s * 2);
  glBegin(GL_POINTS);
  glColor3f(0, 0, 0);
  glVertex3f(x1 * 1.2, y, z1 * 1.2);
  glVertex3f(x2 * 1.2, y, z2 * 1.2);
  glEnd();
  glPopMatrix();
}

// draws a long-haired furball
// pretty much more CPU intensive than it seems
// thus ultimate furballs are disabled in LOFI mode
void
draw_furball_ultimate(FURBALL* f)
{
  int c, x, y, z, v;
  float *vertex, *basevertex, dir[3], dist, lsize, deviation, *colour, *basecolour, colour_deviation;
  deviation = 1;
  colour_deviation = .2;

  // fetches base vertex buffer for augmentation
  vertex = f->mine->vtx;
  basevertex = f->buf->vtx;
  colour = f->mine->clr;
  basecolour = f->buf->clr;

  // augments the buffer
  for (x = 0; x < f->density; x++)
    {
      for (y = 0; y < f->density; y++)
        {
          for (z = 0; z < f->density; z++)
            {
              for (c = 0; c < TRAIL; c++)
                {
                  dir[0] = (f->trail[c][0] - f->x);
                  dir[1] = (f->trail[c][1] - f->y);
                  dir[2] = (f->trail[c][2] - f->z);
                  for (v = 0; v < 3; v++)
                    vertex[v] = basevertex[v] + (dir[v]);

                  dir[0] = f->trail[c + 1][0] - f->x;
                  dir[1] = f->trail[c + 1][1] - f->y;
                  dir[2] = f->trail[c + 1][2] - f->z;
                  for (v = 0; v < 3; v++)
                    vertex[3 + v] = basevertex[3 + v] + (dir[v]);
                  vertex += 6;
                  basevertex += 6;
                  colour += 6;
                  basecolour += 6;
                }
            }
        }
    }

  // calculates line width and draws

  dist = sqrt((f->x - playerx) * (f->x - playerx) + (f->y - playery) * (f->y - playery) + (f->z - playerz) * (f->z - playerz));
  lsize = 640 / dist;
  // if it's too small, no use drawing
  if (lsize < .01)
    return;
  glLineWidth(lsize);
  draw_buffer(f->mine, f->x, f->y, f->z, f->scale, 0);
}

// draws a simple furball from a vertex buffer
// stretches it a bit according to vertical velocity
// and displaces a bit (actually eyes look displaced this way)
void
draw_furball_normal(FURBALL* f)
{
  float stretch, displace;

  stretch = 1.0 + pow(ABS(f->y - f->trail[1][1]) * .1, 2);
  displace = f->y - f->trail[1][1];
  if (stretch > 3.0f)
    stretch = 3.0f;

  draw_buffer_ex(f->mine, f->x, f->y - displace, f->z, f->scale, f->scale * stretch, f->scale, 1.0f, 0.0f);
}

// draws a blood particles array and makes coffee
void
draw_blood(BLOOD* b, int num)
{
  int c;
  glPointSize(20.0f);
  glBegin(GL_POINTS);
  for (c = 0; c < num; c++)
    {
      glColor3fv(&(b[c].r));
      glVertex3fv(&(b[c].x));
    }
  glEnd();
}

// draws a furball (picks one of the above and checks distance)
void
draw_furball(FURBALL* f)
{
  float x, y, z;
  x = f->x - playerx;
  y = f->y - playery;
  z = f->z - playerz;
  if ((x * x + y * y + z * z) < DRAW_DIST * DRAW_DIST)
    {
      if (f->exists)
        {
          if (f->ultimate)
            draw_furball_ultimate(f);
          else
            draw_furball_normal(f);
          draw_eyes(f);
        }
      if (f->dying)
        {
          draw_blood(f->blood, PARTICLES);
        }
    }
}

// draws all furballs
void
draw_ballz()
{
  int c;
  for (c = 0; c < BALLZ; c++)
    {
      draw_furball(ballz[c]);
    }
  glLineWidth(LINE_WIDTH);
}

// draws all entities
// actually draws only the ones that are close enough
void
draw_ents()
{
  int c;
  float x, z, dist;
  for (c = 0; c < num_ents; c++)
    {
      x = ents[c].x - playerx;
      z = (ents[c].z) - playerz;
      dist = (x * x + z * z) - 32 * 32 * ents[c].s * ents[c].s;
      if (dist < DRAW_DIST * DRAW_DIST)
        {
          draw_buffer(ents[c].buf, ents[c].x, ents[c].y, ents[c].z, ents[c].s, ents[c].a);
        }
    }
}

// apparently, this draws ground.
// ground changes colour dynamically, and is not a textured quad
// but a series of quads, drawn only near the player
// this saves cpu and gpu from dynamic texture updates
// it also draws everything else :P
void
draw_tree()
{
  clock_t clk;
  int c, x, y, pix;
  float fx, fy, fs, r, g, b, xx, yy;
  clk = clock();

  // draws array of ground quads
  glBegin(GL_QUADS);
  for (x = 0; x < ground_bmp->w; x++)
    {
      for (y = 0; y < ground_bmp->h; y++)
        {
          fx = (x * WORLD_SIZE) / ground_bmp->w;
          fy = (y * WORLD_SIZE) / ground_bmp->w;
          fs = WORLD_SIZE / ground_bmp->w;
          xx = fx - playerx;
          yy = fy - playerz;

          // if close enough, draws
          if (xx * xx + yy * yy < DRAW_DIST * DRAW_DIST)
            {
              c = (x + y) & 1;
              pix = getpixel(ground_bmp, x, y);
              r = getr(pix) / 256.0;
              g = getg(pix) / 256.0;
              b = getb(pix) / 256.0;
              if (c)
                {
                  r *= .8;
                  g *= .8;
                  b *= .8;
                }

              glColor3f(r, g, b);
              glVertex3f(fx, 0.0, fy);
              glVertex3f(fx, 0.0, fy + fs);
              glVertex3f(fx + fs, 0.0, fy + fs);
              glVertex3f(fx + fs, 0.0, fy);
            }
        }
    }
  glEnd();
  glMatrixMode(GL_MODELVIEW);

  // draws map object
  draw_ents();

  // draws grass
  draw_buffer(grass, 0, 0, 0, 1, 0);

  // draws furballs
  draw_ballz();

  // and particles
  draw_blood(blood, PARTICLES);
}

////////////////////////////////////////////////////
//////////////////////////////////////////////////// GENERATOR FUNCTIONS
////////////////////////////////////////////////////

// generates a mesh from one bitmap
// by 'spinning it' around Y axis
BUFFER*
generate_cloud_single(const char* filename, int density)
{
  BUFFER* b;
  BITMAP* bmp;
  int x, y, z, d, size_x, size_y, size_z, point_counter, img_x, img_y, col, c_x, c_y, c_z;
  float deviation, *vertex, *colour, colour_deviation;
  printf("Creating cloud from %s...\n", filename);
  deviation = .4f;
  colour_deviation = .03f;
  bmp = load_bmp(filename, 0);
  size_x = bmp->w;
  size_y = bmp->h;
  size_z = bmp->w;
  printf("Dimensions are: %i x %i x %i\n", size_x, size_y, size_z);

  // allocates a too large buffer, or too small, if of higher density
  b = malloc(sizeof(BUFFER));
  b->vtx = malloc(size_x * size_y * size_z * 3 * sizeof(float));
  b->clr = malloc(size_x * size_y * size_z * 3 * sizeof(float));
  b->hardware = 0;
  point_counter = 0;
  vertex = b->vtx;
  colour = b->clr;
  c_x = size_x / 2;
  c_z = size_z / 2;
  c_y = size_y / 2;
  printf("Beggining clouding...\n");
  for (x = 0; x < size_x; x++)
    {
      for (y = 0; y < size_y; y++)
        {
          for (z = 0; z < size_z; z++)
            {
              for (d = 0; d < density; d++)
                {
                  // tests a voxel against image, where y is y
                  // and x is distance from centre
                  // if ok, sets colour and displaces a bit
                  img_x = floor(sqrt((x - c_x) * (x - c_x) + (z - c_z) * (z - c_z))) + c_x;
                  img_y = bmp->h - y - 1;
                  col = getpixel(bmp, img_x, img_y);
                  if (col != makecol(255, 0, 255) && img_x < bmp->w)
                    {
                      vertex[0] = (1.0f * x) + ((((rand() % 512) / 256.0f) - 1.0f) * deviation) - c_x;
                      vertex[1] = (1.0f * y) + ((((rand() % 512) / 256.0f) - 1.0f) * deviation);
                      vertex[2] = (1.0f * z) + ((((rand() % 512) / 256.0f) - 1.0f) * deviation) - c_z;
                      colour[0] = (getr(col) / 256.0f) + ((((rand() % 512) / 256.0f) - 1.0f) * colour_deviation);
                      colour[1] = (getg(col) / 256.0f) + ((((rand() % 512) / 256.0f) - 1.0f) * colour_deviation);
                      colour[2] = (getb(col) / 256.0f) + ((((rand() % 512) / 256.0f) - 1.0f) * colour_deviation);
                      vertex += 3;
                      colour += 3;
                      point_counter++;
                    }
                }
            }
        }
    }
  b->size = point_counter;
  b->mode = GL_POINTS;
  vboize(b);
  printf("Created cloud of %i points!\n", point_counter);
  destroy_bitmap(bmp);
  return b;
}

// generates a mesh from two images
// these are axis projections from front and side
BUFFER*
generate_cloud_double(const char* filename_x, const char* filename_z, int density)
{
  BUFFER* b;
  BITMAP *bmpx, *bmpz;
  int x, y, z, d, size_x, size_y, size_z, point_counter, img_x, img_y;
  int img_z, col1, col2, c_x, c_y, c_z;
  float deviation, *vertex, *colour, colour_deviation;
  deviation = .4f;
  colour_deviation = .03f;
  bmpx = load_bmp(filename_x, 0);
  bmpz = load_bmp(filename_z, 0);
  size_x = bmpx->w;
  size_y = bmpx->h;
  size_z = bmpz->w;

  // another bad malloc
  b = malloc(sizeof(BUFFER));
  b->vtx = malloc(size_x * size_y * size_z * 3 * sizeof(float));
  b->clr = malloc(size_x * size_y * size_z * 3 * sizeof(float));
  b->hardware = 0;
  point_counter = 0;
  vertex = b->vtx;
  colour = b->clr;
  c_x = size_x / 2;
  c_z = size_z / 2;
  c_y = size_y / 2;
  for (x = 0; x < size_x; x++)
    {
      for (y = 0; y < size_y; y++)
        {
          for (z = 0; z < size_z; z++)
            {
              for (d = 0; d < density; d++)
                {
                  // this tests the pixel against two projection images
                  // and sets the colour by averaging from those
                  img_x = x;
                  img_z = z;
                  img_y = bmpx->h - y - 1;
                  col1 = getpixel(bmpx, img_x, img_y);
                  col2 = getpixel(bmpz, img_z, img_y);
                  if (col1 != makecol(255, 0, 255) && col2 != makecol(255, 0, 255) && img_x < bmpx->w && img_z < bmpz->w)
                    {
                      vertex[0] = (1.0f * x) + ((((rand() % 512) / 256.0f) - 1.0f) * deviation) - c_x;
                      vertex[1] = (1.0f * y) + ((((rand() % 512) / 256.0f) - 1.0f) * deviation);
                      vertex[2] = (1.0f * z) + ((((rand() % 512) / 256.0f) - 1.0f) * deviation) - c_z;
                      colour[0] = ((getr(col1) + getr(col2)) / 512.0f) + ((((rand() % 512) / 256.0f) - 1.0f) * colour_deviation);
                      colour[1] = ((getg(col1) + getg(col2)) / 512.0f) + ((((rand() % 512) / 256.0f) - 1.0f) * colour_deviation);
                      colour[2] = ((getb(col1) + getb(col2)) / 512.0f) + ((((rand() % 512) / 256.0f) - 1.0f) * colour_deviation);
                      vertex += 3;
                      colour += 3;
                      point_counter++;
                    }
                }
            }
        }
    }
  b->mode = GL_POINTS;
  b->size = point_counter;
  vboize(b);
  destroy_bitmap(bmpx);
  destroy_bitmap(bmpz);
  printf("Created cloud of %i points!\n", point_counter);
  return b;
}

// generates a mesh from 3 projection images (front, side, top)
BUFFER*
generate_cloud_triple(const char* filename_x, const char* filename_z, const char* filename_y, int density)
{
  BUFFER* b;
  BITMAP *bmpx, *bmpz, *bmpy;
  int x, y, z, d, size_x, size_y, size_z, point_counter, img_x, img_y;
  int img_z, col1, col2, col3, c_x, c_y, c_z;
  float deviation, *vertex, *colour, colour_deviation;
  deviation = .4f;
  colour_deviation = .03f;
  bmpx = load_bmp(filename_x, 0);
  bmpz = load_bmp(filename_z, 0);
  bmpy = load_bmp(filename_y, 0);
  size_x = bmpx->w;
  size_y = bmpx->h;
  size_z = bmpz->w;
  printf("Dimensions are: %i x %i x %i\n", size_x, size_y, size_z);

  // bad malloc Mk.3
  b = malloc(sizeof(BUFFER));
  b->vtx = malloc(size_x * size_y * size_z * 3 * sizeof(float));
  b->clr = malloc(size_x * size_y * size_z * 3 * sizeof(float));
  b->hardware = 0;
  point_counter = 0;
  vertex = b->vtx;
  colour = b->clr;
  c_x = size_x / 2;
  c_z = size_z / 2;
  c_y = size_y / 2;
  for (x = 0; x < size_x; x++)
    {
      for (y = 0; y < size_y; y++)
        {
          for (z = 0; z < size_z; z++)
            {
              for (d = 0; d < density; d++)
                {
                  // similar to the above, tests against 3 images and
                  // averages the colour
                  img_x = x;
                  img_z = z;
                  img_y = bmpx->h - y - 1;
                  col1 = getpixel(bmpx, img_x, img_y);
                  col2 = getpixel(bmpz, img_z, img_y);
                  col3 = getpixel(bmpy, img_x, img_z);
                  if (col3 != makecol(255, 0, 255) && col1 != makecol(255, 0, 255) && col2 != makecol(255, 0, 255) && img_x < bmpx->w && img_z < bmpz->w)
                    {
                      vertex[0] = (1.0f * x) + ((((rand() % 512) / 256.0f) - 1.0f) * deviation) - c_x;
                      vertex[1] = (1.0f * y) + ((((rand() % 512) / 256.0f) - 1.0f) * deviation);
                      vertex[2] = (1.0f * z) + ((((rand() % 512) / 256.0f) - 1.0f) * deviation) - c_z;
                      colour[0] = ((getr(col1) + getr(col2)) / 512.0f) + ((((rand() % 512) / 256.0f) - 1.0f) * colour_deviation);
                      colour[1] = ((getg(col1) + getg(col2)) / 512.0f) + ((((rand() % 512) / 256.0f) - 1.0f) * colour_deviation);
                      colour[2] = ((getb(col1) + getb(col2)) / 512.0f) + ((((rand() % 512) / 256.0f) - 1.0f) * colour_deviation);
                      vertex += 3;
                      colour += 3;
                      point_counter++;
                    }
                }
            }
        }
    }
  b->mode = GL_POINTS;
  b->size = point_counter;
  vboize(b);
  destroy_bitmap(bmpx);
  destroy_bitmap(bmpz);
  printf("Created cloud of %i points!\n", point_counter);

  return b;
}

//  this creates a furball instance
FURBALL*
spawn_furball(float x, float y, float z, BUFFER* b, int density, int ultimate, float red, float green, float blue)
{
  FURBALL* f;
  int bufsize, c;
  float *colour, *basecolour;

  // number of mesh vertices is specified through density
  // it represents, slices, sectors and depth
  bufsize = density * density * density * TRAIL * (ultimate ? 12 : 3) * sizeof(float);

  // allocates memory
  f = malloc(sizeof(FURBALL));
  f->buf = b;
  f->mine = malloc(sizeof(BUFFER));
  f->mine->vtx = malloc(bufsize);
  f->mine->clr = malloc(bufsize);
  f->mine->size = b->size;
  f->mine->hardware = 0;
  f->mine->mode = b->mode;

  // copies the buffer from base for dynamic updates
  memcpy(f->mine->vtx, f->buf->vtx, bufsize);
  memcpy(f->mine->clr, f->buf->clr, bufsize);

  // sets colour
  colour = f->mine->clr;
  basecolour = f->buf->clr;
  for (c = 0; c < f->mine->size; c++)
    {
      colour[0] = basecolour[0] * red;
      colour[1] = basecolour[1] * green;
      colour[2] = basecolour[2] * blue;
      colour += 3;
      basecolour += 3;
    }

  // sets position and other properties
  f->x = x;
  f->y = y;
  f->z = z;
  f->r = red;
  f->g = green;
  f->b = blue;
  f->density = density;
  f->scale = 3;
  f->ultimate = ultimate;
  f->exists = 1;
  return f;
}

// generates a base buffer for hairy furball
BUFFER*
generate_furball_ultimate(float size, int cuts, int density)
{
  BUFFER* b;
  int x, y, z, c, point_counter;
  float fx, fy, fz, *vertex, *colour;
  float deviation, colour_deviation;
  deviation = size * .1f;
  colour_deviation = .03;
  float hair_colour, hair_distance;

  // allocates memory
  b = malloc(sizeof(BUFFER));
  b->vtx = malloc(density * density * density * cuts * 12 * sizeof(float));
  b->clr = malloc(density * density * density * cuts * 12 * sizeof(float));
  b->hardware = 0;
  point_counter = 0;
  vertex = b->vtx;
  colour = b->clr;
  for (x = 0; x < density; x++)
    {
      for (y = 0; y < density; y++)
        {
          for (z = 0; z < density; z++)
            {
              // generates a sphere of vertices
              fx = (x * M_PI * 2) / (density - 1);
              fy = (y * M_PI) / (density - 1);
              fz = (z * M_PI * 2) / (density - 1);

              // generates protruding hair
              for (c = 0; c < cuts; c++)
                {

                  hair_colour = (.5 + .5 * (c * 1) / (cuts));
                  hair_distance = size + (c * size) / (cuts);

                  // generates position only for the hair beginning
                  // rest is determined by velocity
                  // this generates startpoint
                  if (!c)
                    {
                      vertex[0] = hair_distance * cos(fx) * sin(fy) + ((((rand() % 512) / 256.0f) - 1.0f) * deviation);
                      vertex[1] = hair_distance * cos(fy) + ((((rand() % 512) / 256.0f) - 1.0f) * deviation);
                      vertex[2] = hair_distance * sin(fz) * sin(fy) + ((((rand() % 512) / 256.0f) - 1.0f) * deviation);
                      colour[0] = hair_colour + ((((rand() % 512) / 256.0f) - 1.0f) * colour_deviation);
                      colour[1] = hair_colour + ((((rand() % 512) / 256.0f) - 1.0f) * colour_deviation);
                      colour[2] = hair_colour + ((((rand() % 512) / 256.0f) - 1.0f) * colour_deviation);
                    }
                  else
                    {
                      vertex[0] = *(vertex - 3);
                      vertex[1] = *(vertex - 2);
                      vertex[2] = *(vertex - 1);
                      colour[0] = *(colour - 3);
                      colour[1] = *(colour - 2);
                      colour[2] = *(colour - 1);
                    }

                  // sets up hair endpoint
                  // hair is actually line list as opposed to a strip
                  hair_colour = ((c + 1) * 1.0f) / (cuts);
                  hair_distance = size + ((c + 1) * size) / (cuts);
                  vertex[3] = hair_distance * cos(fx) * sin(fy) + ((((rand() % 512) / 256.0f) - 1.0f) * deviation);
                  vertex[4] = hair_distance * cos(fy) + ((((rand() % 512) / 256.0f) - 1.0f) * deviation);
                  vertex[5] = hair_distance * sin(fz) * sin(fy) + ((((rand() % 512) / 256.0f) - 1.0f) * deviation);
                  colour[3] = hair_colour + ((((rand() % 512) / 256.0f) - 1.0f) * colour_deviation);
                  colour[4] = hair_colour + ((((rand() % 512) / 256.0f) - 1.0f) * colour_deviation);
                  colour[5] = hair_colour + ((((rand() % 512) / 256.0f) - 1.0f) * colour_deviation);

                  // witty pointer iterators
                  vertex += 6;
                  colour += 6;
                  point_counter += 2;
                }
            }
        }
    }
  b->mode = GL_LINES;
  b->size = point_counter;
  printf("Created furball of %i points!\n", point_counter);

  return b;
}

// egnerates a simple furball
BUFFER*
generate_furball_normal(float size, int cuts, int density)
{
  BUFFER* b;
  int x, y, z, c, point_counter;
  float fx, fy, fz, *vertex, *colour;
  float deviation, colour_deviation, colour_factor;
  deviation = size * .3f;
  colour_deviation = .08;

  // allocates buffers
  b = malloc(sizeof(BUFFER));
  b->vtx = malloc(density * density * density * cuts * 3 * sizeof(float));
  b->clr = malloc(density * density * density * cuts * 3 * sizeof(float));
  b->hardware = 0;
  point_counter = 0;
  vertex = b->vtx;
  colour = b->clr;
  for (x = 0; x < density; x++)
    {
      for (y = 0; y < density; y++)
        {
          for (z = 0; z < density; z++)
            {
              // gnerates vertex position in a sphere, and displaces a bit
              fx = (x * M_PI * 2) / (density - 1);
              fy = (y * M_PI) / (density - 1);
              fz = (z * M_PI * 2) / (density - 1);
              for (c = 0; c < cuts; c++)
                {

                  vertex[0] = size * cos(fx) * sin(fy) + ((((rand() % 512) / 256.0f) - 1.0f) * deviation);
                  vertex[1] = size * cos(fy) + ((((rand() % 512) / 256.0f) - 1.0f) * deviation);
                  vertex[2] = size * sin(fz) * sin(fy) + ((((rand() % 512) / 256.0f) - 1.0f) * deviation);

                  colour_factor = sqrt(vertex[0] * vertex[0] + vertex[1] * vertex[1] + vertex[2] * vertex[2]) / size;
                  colour[0] = colour_factor + ((((rand() % 512) / 256.0f) - 1.0f) * colour_deviation);
                  colour[1] = colour_factor + ((((rand() % 512) / 256.0f) - 1.0f) * colour_deviation);
                  colour[2] = colour_factor + ((((rand() % 512) / 256.0f) - 1.0f) * colour_deviation);
                  vertex += 3;
                  colour += 3;
                  point_counter++;
                }
            }
        }
    }
  b->mode = GL_POINTS;
  b->size = point_counter;
  printf("Created furball of %i points!\n", point_counter);

  return b;
}

// creates all the furballs
void
create_ballz()
{
  int c, is_good;
  float x, y, z, r, g, b;
  for (c = 0; c < BALLZ; c++)
    {
      // if is_good, then it's 'ultimate' == hairy
      is_good = !(rand() % 8);

      // no hairy furballs for lofi version
      if (LOFI)
        is_good = 0;

      // randomises colour and position
      r = ((rand() % 256) / 256.0f);
      g = ((rand() % 256) / 256.0f);
      b = ((rand() % 256) / 256.0f);
      x = ((rand() % 1024) / 1024.0f) * WORLD_SIZE;
      y = HEIGHT;
      z = ((rand() % 1024) / 1024.0f) * WORLD_SIZE;

      // creates the instance
      ballz[c] = spawn_furball(x, y, z, is_good ? ultimate_furball : casual_furball, 8, is_good, r, g, b);

      // sets up ai
      ballz[c]->iq = 0;
      ballz[c]->smart = .05f + ((rand() % 256) / 256.0f) * .15;

      // adjusts ai for hairy furball
      if (is_good)
        ballz[c]->smart *= .2;
      if (!is_good)
        ballz[c]->scale = .2f + ((rand() % 256) / 256.0f) * .15;
      else
        ballz[c]->scale = 1.0f + ((rand() % 256) / 256.0f) * 3.0f;
      ballz[c]->dying = 0;
    }
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////  MECHANINCS FUNCTIONS
////////////////////////////////////////////////////

// creates [num] blood particles at [b] buffer
void
explode(BLOOD* b, int num, float x, float y, float z)
{
  int c;
  for (c = 0; c < num; c++)
    {
      // sets position
      b[c].x = x;
      b[c].y = y;
      b[c].z = z;
      // random velocity and colour
      b[c].vx = ((((rand() % 512) / 256.0f) - 1.0f) * 4.0f);
      b[c].vy = ((((rand() % 512) / 256.0f) - 1.0f) * 8.0f);
      b[c].vz = ((((rand() % 512) / 256.0f) - 1.0f) * 4.0f);
      b[c].r = .6 + ((((rand() % 512) / 256.0f) - 1.0f) * .4f);
      b[c].g = 0;
      b[c].b = 0;
      b[c].alive = 1;
    }
}

// shoots gun
int
shoot()
{
  int c, h = 0;
  FURBALL* hitball = 0;
  float p[3], b[3], m[3], p_minus_b[3], t0, hit[3], pdist, rdist, rhit[3];
  play_sample(shot, 255, 128, 1000, 0);
  // gets shoot position and direction
  b[0] = playerx;
  b[1] = playery;
  b[2] = playerz;
  m[0] = cos(looky) * cos(lookx);
  m[1] = -sin(lookx);
  m[2] = sin(looky) * cos(lookx);
  pdist = 999999.0f;
  for (c = 0; c < BALLZ; c++)
    {
      // picks furball closest to the ray
      if (ballz[c]->exists)
        {
          p[0] = ballz[c]->x;
          p[1] = ballz[c]->y;
          p[2] = ballz[c]->z;
          p_minus_b[0] = p[0] - b[0];
          p_minus_b[1] = p[1] - b[1];
          p_minus_b[2] = p[2] - b[2];
          t0 = dot3v(m, p_minus_b) / dot3v(m, m);
          if (t0 > 0)
            {
              rhit[0] = p[0] - (b[0] + t0 * m[0]);
              rhit[1] = p[1] - (b[1] + t0 * m[1]);
              rhit[2] = p[2] - (b[2] + t0 * m[2]);
              rdist = length3v(rhit);
              if (rdist < 10.0f)
                {
                  printf("hit with rdist: %f\n", rdist);
                  if (dist3v(b, p) < pdist)
                    {
                      pdist = dist3v(b, p);
                      hit[0] = p[0];
                      hit[1] = p[1];
                      hit[2] = p[2];
                      h = 1;
                      hitball = ballz[c];
                    }
                }
            }
        }
    }
  // checks if furball is close enough to the ray to be hit
  if (h && pdist < 300.0f)
    {
      // disables furball rendering
      hitball->exists = 0;

      // enables blood renderign for 100 frames
      hitball->dying = 100;

      // creates blood
      explode(hitball->blood, PARTICLES, hitball->x, hitball->y, hitball->z);

      // plays sound
      play_sample(slash, 255, 128, 900 + (rand() % 200), 0);
      play_sample(die[rand() % 8], 255, 128, 1000, 0);
    }
  return h;
}

// updates blood particles
void
update_blood(BLOOD* b, int num)
{
  int c, pic_x, pic_y, pr, pg, pb, pix;
  for (c = 0; c < num; c++)
    {
      // if the particle is alive...
      if (b[c].alive)
        {

          // updates position
          b[c].x += b[c].vx;
          b[c].y += b[c].vy;
          b[c].z += b[c].vz;

          // adds gravity
          b[c].vy -= .6;

          // if it hits the ground
          if (b[c].y < 0)
            {
              // no more updates for this one
              b[c].alive = 0;

              // if was inside world sware, colours a pixel red
              pic_x = (b[c].x * ground_bmp->w) / (WORLD_SIZE);
              pic_y = (b[c].z * ground_bmp->h) / (WORLD_SIZE);
              if (pic_x >= 0 && pic_y >= 0 && pic_x < ground_bmp->w && pic_y < ground_bmp->h)
                {
                  pix = getpixel(ground_bmp, pic_x, pic_y);
                  pr = (getr(pix) + b[c].r * 128) / 1.5;
                  pg = (getg(pix) + b[c].g * 128) / 1.5;
                  pb = (getb(pix) + b[c].b * 128) / 1.5;
                  putpixel(ground_bmp, pic_x, pic_y, makecol(pr, pg, pb));
                }
            }
        }
    }
}

// updates furballs
void
update_ballz()
{
  int c, i;
  FURBALL* fur;
  for (c = 0; c < BALLZ; c++)
    {
      // if it's still alive
      if (ballz[c]->exists)
        {
          fur = ballz[c];

          // decreases direction change counter
          // and changes movement direction if needed
          fur->iq -= fur->smart;
          if (fur->iq <= 0)
            {
              fur->iq = (IQ + (rand() % IQ)) * M_PI;
              fur->speed = ((rand() % 256) / 256.0f) * SPEED;
              fur->a = ((rand() % 256) / 256.0f) * M_PI * 2;
              fur->bounce_rate = 10.0f + ((rand() % 256) / 256.0f) * 40.0f;
            }

          // adjusts bounce rate
          fur->bounce += .1f * SPEED;

          // updates position
          fur->x += cos(fur->a) * fur->speed;
          fur->z += sin(fur->a) * fur->speed;
          fur->y = fur->scale * (fur->ultimate ? 1 : 8) + ABS(sin(fur->iq) * fur->bounce_rate);

          // keeps inside the world
          if (fur->x < 0)
            fur->x = 0;
          if (fur->z < 0)
            fur->z = 0;
          if (fur->x > WORLD_SIZE)
            fur->x = WORLD_SIZE;
          if (fur->z > WORLD_SIZE)
            fur->z = WORLD_SIZE;

          // wites down trail for hair rendering
          for (i = TRAIL; i > 0; i--)
            {
              fur->trail[i][0] = fur->trail[i - 1][0];
              fur->trail[i][1] = fur->trail[i - 1][1];
              fur->trail[i][2] = fur->trail[i - 1][2];
            }
          fur->trail[0][0] = fur->x;
          fur->trail[0][1] = fur->y;
          fur->trail[0][2] = fur->z;
        }
      else if (ballz[c]->dying)
        {
          // if it's dying, only particles get updated
          update_blood(ballz[c]->blood, PARTICLES);
          ballz[c]->dying--;
        }
    }
}

// generates world map from 3 files
// first is colour, where it stores ground and grass colours
// second is grass height
// third is entity distribution, a specific colour represents one entity
// made so for quick level editing
BUFFER*
generate_world_map(const char* grass_file, const char* height_file, const char* entity_file, int density)
{
  BUFFER* b;
  BITMAP *bmpg, *bmph, *bmph_t, *bmpg_t;
  int x, y, d, size_x, size_y, line_counter, col;
  float deviation, *vertex, *colour, colour_deviation, grass_x, grass_y, grass_h, patch_size;
  printf("Creating grass from %s...\n", grass_file);

  // randomisation constants
  deviation = .4f;
  colour_deviation = .1f;

  // lofi gets smaller bitmaps
  if (LOFI)
    {
      bmpg_t = load_bmp(grass_file, 0);
      bmph_t = load_bmp(height_file, 0);
      bmpg = create_bitmap(bmpg_t->w / 2, bmpg_t->h / 2);
      bmph = create_bitmap(bmph_t->w / 2, bmph_t->h / 2);
      stretch_blit(bmph_t, bmph, 0, 0, bmph_t->w, bmph_t->h, 0, 0, bmph_t->w / 2, bmph_t->h / 2);
      stretch_blit(bmpg_t, bmpg, 0, 0, bmpg_t->w, bmpg_t->h, 0, 0, bmpg_t->w / 2, bmpg_t->h / 2);
      destroy_bitmap(bmph_t);
      destroy_bitmap(bmpg_t);
    }
  else
    {
      bmpg = load_bmp(grass_file, 0);
      bmph = load_bmp(height_file, 0);
    }
  size_x = bmpg->w;
  size_y = bmpg->h;
  printf("Dimensions are: %i x %i x %i\n", size_x, size_y, density);

  // allocates grass buffer
  b = malloc(sizeof(BUFFER));
  b->vtx = malloc(size_x * size_y * density * 6 * sizeof(float));
  b->clr = malloc(size_x * size_y * density * 6 * sizeof(float));
  b->hardware = 0;
  line_counter = 0;
  vertex = b->vtx;
  colour = b->clr;
  patch_size = WORLD_SIZE / bmpg->w;
  printf("Beggining clouding...\n"); // actually grassing
  for (x = 0; x < size_x; x++)
    {
      for (y = 0; y < size_y; y++)
        {
          for (d = 0; d < density; d++)
            {
              // for each bitmap pixel it creates [density] grass blades
              // colour and height get randomised a bit
              grass_x = ((x * WORLD_SIZE) / size_x) + ((((rand() % 512) / 512.0f)) * patch_size);
              grass_y = ((y * WORLD_SIZE) / size_x) + ((((rand() % 512) / 512.0f)) * patch_size);
              grass_h = (.1 * getr(getpixel(bmph, x, y))) + ((((rand() % 512) / 512.0f)) * (grass_h * .3f));
              vertex[0] = grass_x;
              vertex[1] = 0;
              vertex[2] = grass_y;
              vertex[3] = grass_x;
              vertex[4] = grass_h;
              vertex[5] = grass_y;
              col = getpixel(bmpg, x, y);
              colour[3] = ((getr(col)) / 256.0f) + ((((rand() % 512) / 256.0f) - 1.0f) * colour_deviation);
              colour[4] = ((getg(col)) / 256.0f) + ((((rand() % 512) / 256.0f) - 1.0f) * colour_deviation);
              colour[5] = ((getb(col)) / 256.0f) + ((((rand() % 512) / 256.0f) - 1.0f) * colour_deviation);
              colour[0] = colour[3] * .5f;
              colour[1] = colour[4] * .5f;
              colour[2] = colour[5] * .5f;
              line_counter += 2;
              vertex += 6;
              colour += 6;
            }
        }
    }
  printf("Generated %i lines!\n", line_counter);
  b->mode = GL_LINES;
  b->size = line_counter;
  vboize(b);
  ground_bmp = bmpg;

  destroy_bitmap(bmph);

  // creates entities
  bmph = load_bmp(entity_file, 0);
  size_x = bmph->w;
  size_y = bmph->h;
  for (x = 0; x < size_x; x++)
    {
      for (y = 0; y < size_y; y++)
        {
          // iterates through pixel looking for values representing entities
          grass_x = ((x * WORLD_SIZE) / size_x) + ((((rand() % 512) / 512.0f)) * patch_size);
          grass_y = ((y * WORLD_SIZE) / size_x) + ((((rand() % 512) / 512.0f)) * patch_size);
          grass_h = 1.0f + ((((rand() % 512) / 512.0f)) * 1.0f);
          col = getpixel(bmph, x, y);
          if (col == makecol(0, 255, 0))
            {
              ents[num_ents].buf = tall_tree;
              ents[num_ents].x = grass_x;
              ents[num_ents].y = 0;
              ents[num_ents].z = grass_y;
              ents[num_ents].s = grass_h * 4;
              ents[num_ents].a = ((((rand() % 512) / 512.0f)) * 360);
              num_ents++;
            }
          if (col == makecol(255, 255, 128))
            {
              ents[num_ents].buf = cactus;
              ents[num_ents].x = grass_x;
              ents[num_ents].y = 0;
              ents[num_ents].z = grass_y;
              ents[num_ents].s = grass_h * 2;
              ents[num_ents].a = ((((rand() % 512) / 512.0f)) * 360);
              num_ents++;
            }
          if (col == makecol(0, 128, 0))
            {
              ents[num_ents].buf = palm;
              ents[num_ents].x = grass_x;
              ents[num_ents].y = 0;
              ents[num_ents].z = grass_y;
              ents[num_ents].s = grass_h * 3;
              ents[num_ents].a = ((((rand() % 512) / 512.0f)) * 360);
              num_ents++;
            }
          if (col == makecol(128, 64, 64))
            {
              ents[num_ents].buf = stick;
              ents[num_ents].x = grass_x;
              ents[num_ents].y = 0;
              ents[num_ents].z = grass_y;
              ents[num_ents].s = grass_h * 1.5;
              ents[num_ents].a = ((((rand() % 512) / 512.0f)) * 360);
              num_ents++;
            }
          if (col == makecol(0, 255, 255))
            {
              ents[num_ents].buf = pine;
              ents[num_ents].x = grass_x;
              ents[num_ents].y = 0;
              ents[num_ents].z = grass_y;
              ents[num_ents].s = grass_h * 3;
              ents[num_ents].a = ((((rand() % 512) / 512.0f)) * 360);
              num_ents++;
            }
          if (col == makecol(255, 0, 0))
            {
              ents[num_ents].buf = hut;
              ents[num_ents].x = grass_x;
              ents[num_ents].y = 0;
              ents[num_ents].z = grass_y;
              ents[num_ents].s = grass_h * 2;
              ents[num_ents].a = (((rand() % 4)) * 90.0f);
              num_ents++;
            }
          if (col == makecol(255, 128, 128))
            {
              ents[num_ents].buf = church;
              ents[num_ents].x = grass_x;
              ents[num_ents].y = 0;
              ents[num_ents].z = grass_y;
              ents[num_ents].s = grass_h * 6;
              ents[num_ents].a = 90.0f; //((((rand()%512)/512.0f))*360);
              num_ents++;
            }
          if (col == makecol(128, 0, 0))
            {
              ents[num_ents].buf = brickhouse;
              ents[num_ents].x = grass_x;
              ents[num_ents].y = 0;
              ents[num_ents].z = grass_y;
              ents[num_ents].s = grass_h * 1.5;
              ents[num_ents].a = (((rand() % 4)) * 90.0f);
              num_ents++;
            }
        }
    }
  destroy_bitmap(bmph);
  return b;
}

// generate everything batch
void
generate_stuff()
{
  // generate_font("cartographer.ttf",1000);
  ultimate_furball = generate_furball_ultimate(1.0f, TRAIL, 8);
  casual_furball = generate_furball_normal(8.0f, TRAIL, 8);
  tall_tree = generate_cloud_triple("tall_tree.bmp", "tall_tree.bmp", "tall_tree_top.bmp", 1);
  withered_bush = generate_cloud_double("small_bush_left.bmp", "small_bush_front.bmp", 1);
  cactus = generate_cloud_double("cactus_front.bmp", "cactus_side.bmp", 1);
  palm = generate_cloud_single("palm.bmp", 1);
  stick = generate_cloud_double("stick_front.bmp", "stick_side.bmp", 1);
  pine = generate_cloud_single("pine.bmp", 1);
  hut = generate_cloud_double("hut_front.bmp", "hut_side.bmp", 1);
  fence = generate_cloud_triple("fence_side.bmp", "fence_side.bmp", "fence_top.bmp", 3);
  church = generate_cloud_double("church_front.bmp", "church_side.bmp", 1);
  brickhouse = generate_cloud_double("brickhouse_front.bmp", "brickhouse_side.bmp", 1);
  grass = generate_world_map("grassmap.bmp", "grassheight.bmp", "entitymap.bmp", LOFI ? 2 : 4);
  create_ballz();
}

// timed function that updates everything
// and handles input
void
timer_proc(void)
{
  int c;
  float move_spd = MOVE_SPEED;
  clock_t clk;
  clk = clock();

  // if it's intro, waits for space
  if (game_state == INTRO)
    {
      position_mouse(SCREEN_W / 2, SCREEN_H / 2);
      if (key[KEY_SPACE])
        {
          play_sample(youknow, 255, 128, 1000, 0);
          game_state = PLAY;
        }
    }
  else
    {
      // furball talking got disabled
      if (talking)
        {
          talk_counter--;
          if (talk_counter <= 0)
            {
              talk_counter = TALK_DELAY + (rand() % TALK_DELAY);
              // play_sample(furtalk[(rand()%16)>>1],128,128,1000,0);
            }
        }

      // furball whistle makes all furballs
      // move towards player for some time
      if (whistle_timeout)
        {
          whistle_timeout--;
        }
      else
        {
          if (key[KEY_SPACE])
            {
              whistle_timeout = 100;
              play_sample(whistle, 255, 128, 1000, 0);
              for (c = 0; c < BALLZ; c++)
                {
                  ballz[c]->a = atan2(-ballz[c]->z + playerz, -ballz[c]->x + playerx);
                  ballz[c]->iq += M_PI * 4;
                }
            }
        }
      // shift makes you run faster
      if (key[KEY_LSHIFT])
        move_spd *= 2;

      // rotates camera using mouse input
      get_mouse_mickeys(&mx, &my);
      position_mouse(SCREEN_W / 2, mouse_y);
      update_blood(blood, PARTICLES);
      lookx = (M_PI / 2) * (1.0f * (mouse_y - (SCREEN_H / 2))) / (SCREEN_H / 2.0f);
      looky += (M_PI) * (mx * 1.0f) / (640 / 2.0f);

      // moves player according to key input
      if (key[KEY_LEFT] || key[KEY_A])
        {
          playerx += sin(looky) * move_spd;
          playerz += -cos(looky) * move_spd;
          state = WALK;
        }
      if (key[KEY_RIGHT] || key[KEY_D])
        {
          playerx += -sin(looky) * move_spd;
          playerz += cos(looky) * move_spd;
          state = WALK;
        }

      if (key[KEY_UP] || key[KEY_W])
        {
          playerx += cos(looky) * move_spd;
          playerz += sin(looky) * move_spd;
          state = WALK;
        }
      if (key[KEY_DOWN] || key[KEY_S])
        {
          playerx -= cos(looky) * move_spd;
          playerz -= sin(looky) * move_spd;
          state = WALK;
        }

      // if not walking, sets idle state
      if (!key[KEY_UP] && !key[KEY_DOWN] && state != IDLE)
        state = IDLE;

      // if not moving, add breathing effect
      if (state == STOP)
        {
          if (ABS(playery - HEIGHT) > HAYSTACK * sin((playerx / WORLD_SIZE) * M_PI) * sin((playerx / WORLD_SIZE) * M_PI) * sin((playerz / WORLD_SIZE) * M_PI) * sin((playerz / WORLD_SIZE) * M_PI))
            {
              playery -= SGN(playery - HEIGHT) * .4;
            }
          else
            {
              state = IDLE;
              lookf = 0.0;
            }
        }
      // shoots
      if (mouse_b)
        shoot();
    }

  bounce = BOUNCE_RATE;

  // if idle, add breathe effect aswell
  if (state == IDLE)
    {
      lookf += .2;
      playery = HEIGHT + (bounce * .3 * (sin(lookf)));
    }
  if (state == WALK)
    {
      lookf += 0.1;

      playery = HEIGHT + (bounce * (sin(lookf)));
    }

  // keeps player inside the world
  if (playerx > (WORLD_SIZE - 50.0))
    playerx = (WORLD_SIZE - 50.0);
  if (playerz > (WORLD_SIZE - 50.0))
    playerz = (WORLD_SIZE - 50.0);
  if (playerx < 50.0)
    playerx = 50.0;
  if (playerz < 50.0)
    playerz = 50.0;

  // updates furballs
  update_ballz();
}

// rendering function
void
draw(void)
{
  clock_t clk;
  int c;
  float mx, my;
  float ns;
  char nums[16];
  int left = 0;

  // checks for number of furballs left
  for (c = 0; c < BALLZ; c++)
    {
      if (ballz[c]->exists)
        left++;
    }

  // a cheatcode ;)
  if (key[KEY_Y])
    left = 0;

  // sets state to outro if all furballs are dead
  if (!left)
    game_state = OUTRO;

  // stop furball talking if you killed too much
  if (left < 100)
    talking = 0;
  sprintf(nums, "%03i:%03i", left, BALLZ);

  // up sets furball count display
  for (c = 0; c < 7; c++)
    {
      nums[c] -= 0x30;
    }
  clk = clock();

  // basic gl frame setup
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glPushMatrix();
  glFrustum(-((1.0 * SCREEN_W) / SCREEN_H), (1.0 * SCREEN_W) / SCREEN_H, -1.0, 1.0, 1.0, DRAW_DIST);
  glMatrixMode(GL_MODELVIEW);

  // Clear the RGB buffer and the depth buffer
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Set the modelview matrix to be the identity matrix
  glLoadIdentity();
  glPushMatrix();

  // Set the camera
  glRotatef(DEG(lookx), 1, 0, 0);
  glRotatef(DEG(looky) + 90.0, 0, 1, 0);
  glRotatef(DEG(lookz), 0, 0, 1);
  glTranslatef(-playerx, -playery, -playerz);

  // Save the camera matrix
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glTranslatef(0, sin(lookx * 1.333), 0);
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glDisable(GL_DEPTH_TEST);
  glDepthMask(GL_FALSE);

  // draws backgorund
  glBegin(GL_QUADS);
  glColor3f(1.0f, .8, .4);
  glVertex2f(-1, .5);
  glVertex2f(1, .5);
  glColor3f(1.0f, .85, .6);
  glVertex2f(1, .2);
  glVertex2f(-1, .2);

  glVertex2f(-1, .2);
  glVertex2f(1, .2);
  glColor3f(1.0f, .9, .8);
  glVertex2f(1, 0);
  glVertex2f(-1, 0);

  glVertex2f(-1, 0);
  glVertex2f(1, 0);
  glColor3f(0, 0, 0);
  glVertex2f(1, -1);
  glVertex2f(-1, -1);
  glEnd();

  glDepthMask(GL_TRUE);
  glEnable(GL_DEPTH_TEST);
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();

  // Translate and rotate the object
  glColor3f(1.0, 1.0, 1.0);

  // draws world
  draw_tree();

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();

  glDisable(GL_TEXTURE_2D);
  glColor3f(1.0, 1.0, 1.0);
  mx = (mouse_x - 160) / 160.0;
  my = -(mouse_y - 120) / 120.0;

  glPointSize(3.0);
  glDisable(GL_DEPTH_TEST);
  if (game_state == PLAY)
    {

      // draws crosshair
      glBegin(GL_POINTS);
      glColor4f(1.0, 1.0, 1.0, 1.0);
      glVertex2f(0, 0);
      glEnd();

      // draws numbers
      glEnable(GL_TEXTURE_2D);
      glBindTexture(GL_TEXTURE_2D, numbers);
      glBegin(GL_QUADS);

      for (c = 0, ns = 0; c < 7; c++, ns += NUMSIZE)
        {

          glTexCoord2f(nums[c] * .090909f, 1.0);
          glVertex2f(-1.0 + ns, 1.0);
          glTexCoord2f((nums[c] + 1) * .090909f, 1.0);
          glVertex2f(-1.0 + ns + (NUMSIZE * 1.2), 1.0);
          glTexCoord2f((nums[c] + 1) * .090909f, 0.0);
          glVertex2f(-1.0 + ns + (NUMSIZE * 1.2), 1.0 - NUMSIZE);
          glTexCoord2f(nums[c] * .090909f, 0.0);
          glVertex2f(-1.0 + ns, 1.0 - NUMSIZE);
        }

      glEnd();
      glDisable(GL_TEXTURE_2D);
    }

  // renders intro screen
  if (game_state == INTRO)
    {
      glEnable(GL_TEXTURE_2D);
      glBindTexture(GL_TEXTURE_2D, intro);
      glBegin(GL_QUADS);
      glTexCoord2f(0.0, 1.0);
      glVertex2f(-1.0, 1.0);
      glTexCoord2f(1.0, 1.0);
      glVertex2f(1.0, 1.0);
      glTexCoord2f(1.0, 0.0);
      glVertex2f(1.0, -1.0);
      glTexCoord2f(0.0, 0.0);
      glVertex2f(-1.0, -1.0);
      glEnd();
      glDisable(GL_TEXTURE_2D);
    }

  // renders outro screen
  if (game_state == OUTRO)
    {
      glEnable(GL_TEXTURE_2D);
      glBindTexture(GL_TEXTURE_2D, outro);
      glBegin(GL_QUADS);
      glTexCoord2f(0.0, 1.0);
      glVertex2f(-1.0, 1.0);
      glTexCoord2f(1.0, 1.0);
      glVertex2f(1.0, 1.0);
      glTexCoord2f(1.0, 0.0);
      glVertex2f(1.0, -1.0);
      glTexCoord2f(0.0, 0.0);
      glVertex2f(-1.0, -1.0);
      glEnd();
      glDisable(GL_TEXTURE_2D);
    }
  glEnable(GL_DEPTH_TEST);

  // flip buffers
  allegro_gl_flip();
  // printf("Draw took %u msec\n",clock()-clk);
}

// timer thread
void
timer(void)
{
  tim++;
}
END_OF_FUNCTION(timer);

// application entry point
int
main(int argc, char** argv)
{
  int c;
  int w, h;
  RECT rect;
  float fogc[] = {1.0f, .8, .4, .0};
  char samp[64];

  h = rect.bottom - rect.top;
  w = rect.right - rect.left;
  printf("Width: %i\nHeight: %i", w, h);
  // allegro setup
  allegro_init();
  install_allegro_gl();

  allegro_gl_clear_settings();
  allegro_gl_set(AGL_COLOR_DEPTH, 32);
  allegro_gl_set(AGL_Z_DEPTH, 24);
  allegro_gl_set(AGL_WINDOWED, FALSE);
  allegro_gl_set(AGL_DOUBLEBUFFER, 1);
  allegro_gl_set(AGL_SUGGEST, AGL_COLOR_DEPTH | AGL_Z_DEPTH
                                | AGL_DOUBLEBUFFER | AGL_FULLSCREEN);

  if (set_gfx_mode(GFX_OPENGL, w, h, 0, 0) < 0)
    {
      allegro_message("Error setting OpenGL graphics mode:\n%s\n"
                      "Allegro GL error : %s\n",
                      allegro_error, allegro_gl_error);
      return -1;
    }

  install_keyboard();
  install_timer();
  install_mouse();
  install_sound(DIGI_AUTODETECT, MIDI_NONE, argv[0]);

  LOCK_FUNCTION(timer);
  LOCK_VARIABLE(tim);

  // GL setup
  glShadeModel(GL_SMOOTH);
  glEnable(GL_DEPTH_TEST);
  glCullFace(GL_FRONT_AND_BACK);
  glDisable(GL_CULL_FACE);
  glLineWidth(12.0f);
  glPointSize(26.0f);
  glHint(GL_LINE_SMOOTH_HINT, GL_FASTEST);
  glHint(GL_POINT_SMOOTH_HINT, GL_FASTEST);
  glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_FASTEST);
  glHint(GL_POLYGON_SMOOTH_HINT, GL_FASTEST);
  glHint(GL_FOG_HINT, GL_FASTEST);
  glEnable(GL_FOG);
  glFogf(GL_FOG_MODE, GL_LINEAR);

  glFogf(GL_FOG_START, DRAW_DIST * (LOFI ? .8f : .6f));
  glFogf(GL_FOG_END, DRAW_DIST);
  glFogfv(GL_FOG_COLOR, (GLfloat*)&fogc);

  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);

  // timer setup
  install_int_ex(timer, BPS_TO_TIMER(30));
  srand(time(NULL));

  // generates everything
  glEnable(GL_TEXTURE_2D);
  generate_stuff();

  // more gl setup
  glClearColor(1.0f, .8, .4, .0);
  glEnable(GL_BLEND);
  glEnable(GL_ALPHA_TEST);
  glAlphaFunc(GL_GREATER, .3);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  playerx = 20.0;
  playerz = 20.0;
  looky = M_PI / 4;

  get_mouse_mickeys(&mx, &my);

  // loads textures and audio
  intro = SOIL_load_OGL_texture("start.tga", 4, 0, SOIL_FLAG_POWER_OF_TWO | SOIL_FLAG_TEXTURE_REPEATS);
  outro = SOIL_load_OGL_texture("finish.tga", 4, 0, SOIL_FLAG_POWER_OF_TWO | SOIL_FLAG_TEXTURE_REPEATS);
  numbers = SOIL_load_OGL_texture("nums.tga", 4, 0, SOIL_FLAG_POWER_OF_TWO | SOIL_FLAG_TEXTURE_REPEATS);
  youknow = load_wav("youknow.wav");
  slash = load_wav("slash.wav");
  music = load_wav("theme.wav");
  shot = load_wav("shoot.wav");
  whistle = load_wav("whistle.wav");

  for (c = 0; c < SAMPS; c++)
    {
      sprintf(samp, "fur%i.wav", c + 1);
      printf("Loading %s\n", samp);
      furtalk[c] = load_wav(samp);
      sprintf(samp, "die%i.wav", c + 1);
      printf("Loading %s\n", samp);
      die[c] = load_wav(samp);
    }

  // sets player initial position and direction
  playerx = 150;
  playery = HEIGHT;
  playerz = 150;
  lookx = (M_PI / 2) * (1.0f * ((SCREEN_H / 2) - (SCREEN_W / 2))) / (SCREEN_W / 2.0f);
  looky = 1.639519;
  lookz = 0;

  // plays music
  play_sample(music, 255, 128, 1000, 1);

  // main loop runs until someone presses escape
  do
    {
      while (tim > 0)
        {
          timer_proc();
          tim--;
        }
      draw();
      // Sleep(5);
    }
  while (1);//!GetAsyncKeyState(VK_ESCAPE));

  set_gfx_mode(GFX_TEXT, 0, 0, 0, 0);

  return 0;
}
END_OF_MAIN()
