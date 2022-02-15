

#include <allegro.h>
#include <stdio.h>
#include <winalleg.h>
#include <windows.h>

#include "alleggl.h"
#include <SOIL.h>
#include <math.h>

#define LOFI 1
#define TURN_SPEED .1
#define MOVE_SPEED 1.0
#define NUM_TREES 50000
#define WORLD_SIZE 2000.0
#define GRID 256
#define HAYSTACK 100.0
#define FRAND(x) (((rand() % (int)x) * 1024.0) / 1024.0)
#define DEG(n) ((n)*180.0 / M_PI)
//#define ABS(n) ((n)>0 ? (n) : -(n))
//#define SGN(n)  (n>=0 ? 1 : -1)
#define WALK 1
#define IDLE 2
#define STOP 3
#define HEIGHT 10.0
#define HEIGHTF 40.0
#define MAXV 2048
#define NUM_CLOUDS 20
#define CURSOR_SIZE .04
#define WIGGLE 4
#define WORLD_BUFFERS 4096
#define BUFFER_SIZE 300000
#define BOUNCE_RATE 2.0f
#define NUMSIZE .1f

#define FUR_SLICE 8
#define INTRO 1
#define PLAY 2
#define OUTRO 3

#define IQ 2
#define SPEED 1.5f
#if LOFI == 1
#define DRAW_DIST 500.0f
#define PARTICLES 128
#define BALLZ 666
#define TRAIL 1
#define LINE_WIDTH 24.0f
#else
#define LINE_WIDTH 12.0f
#define PARTICLES 256
#define BALLZ 666
#define TRAIL 4
#define DRAW_DIST 2000.0f
#endif
#define MAX_ENTITIES 1024
#define SAMPS 8
#define TALK_DELAY 100

typedef struct BUFFER
{
  float *vtx, *clr;
  GLuint vtx_handle, clr_handle;
  int size, hardware;
  int frames;
  int current_frame;
  GLenum mode;
} BUFFER;

typedef struct BLOOD
{
  float x, y, z;
  float r, g, b, a;
  float vx, vy, vz;
  float s;
  int alive;
} BLOOD;

typedef struct FURBALL
{
  BUFFER *buf, *mine;
  float x, y, z, scale, r, g, b, ultimate, a, bounce, speed, bounce_rate;
  float trail[TRAIL + 1][3];
  int density, exists, dying;
  BLOOD blood[PARTICLES];
  float iq, smart;
} FURBALL;

typedef struct ENTITY
{
  float x, y, z, s, a;
  BUFFER* buf;
} ENTITY;

int frames = 0, state = IDLE;
volatile int tim;
float playerx, playery, playerz, lookx, looky, lookz, lookf;
BUFFER *ultimate_furball, *casual_furball;
BUFFER *grass, *tall_tree, *withered_bush, *cactus, *palm, *stick, *pine, *hut, *fence, *church, *brickhouse;
BITMAP* ground_bmp;
float bounce = BOUNCE_RATE;
FURBALL* ballz[BALLZ];
int mx, my;
BLOOD blood[PARTICLES];
ENTITY ents[MAX_ENTITIES];
int num_ents = 0;
int game_state = INTRO;
GLuint numbers, intro, outro;
SAMPLE *youknow, *furtalk[SAMPS], *die[SAMPS], *slash, *music, *shot, *whistle;
int talk_counter = TALK_DELAY, talking = 1;
int whistle_timeout = 100;

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

float
length3v(float* a)
{
  return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

float
dist3v(float* a, float* b)
{
  float x, y, z;
  x = a[0] - b[0];
  y = a[1] - b[1];
  z = a[2] - b[2];
  return sqrt(x * x + y * y + z * z);
}

float
dot3v(float* a, float* b)
{
  return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

void
explode(BLOOD* b, int num, float x, float y, float z)
{
  int c;
  float ax, ay, az;
  for (c = 0; c < num; c++)
    {
      b[c].x = x;
      b[c].y = y;
      b[c].z = z;
      b[c].vx = ((((rand() % 512) / 256.0f) - 1.0f) * 4.0f);
      b[c].vy = ((((rand() % 512) / 256.0f) - 1.0f) * 8.0f);
      b[c].vz = ((((rand() % 512) / 256.0f) - 1.0f) * 4.0f);
      b[c].r = .6 + ((((rand() % 512) / 256.0f) - 1.0f) * .4f);
      b[c].g = 0;
      b[c].b = 0;
      b[c].alive = 1;
    }
}

int
shoot()
{
  int c, h = 0;
  FURBALL* hitball = 0;
  float p[3], b[3], m[3], p_minus_b[3], t0, hit[3], pdist, rdist, rhit[3];
  play_sample(shot, 255, 128, 1000, 0);
  b[0] = playerx;
  b[1] = playery;
  b[2] = playerz;
  m[0] = cos(looky) * cos(lookx);
  m[1] = -sin(lookx);
  m[2] = sin(looky) * cos(lookx);
  pdist = 999999.0f;
  for (c = 0; c < BALLZ; c++)
    {
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
              if (rdist < 10.0f) //(ballz[c]->scale *(ballz[c]->ultimate ? 1 : 8)))
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
  if (h && pdist < 300.0f)
    {
      printf("Hit at %f %f %f dist: %f\n", hit[0], hit[1], hit[2], pdist);
      // explode(blood,PARTICLES,hit[0],hit[1],hit[2]);
      hitball->exists = 0;
      hitball->dying = 100;
      explode(hitball->blood, PARTICLES, hitball->x, hitball->y, hitball->z);
      play_sample(slash, 255, 128, 900 + (rand() % 200), 0);
      play_sample(die[rand() % 8], 255, 128, 1000, 0);
    }
  else
    printf("No hit!\n");
  return h;
}

void
vboize(BUFFER* b)
{
  // return;
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

void
draw_buffer_ex(BUFFER* b, float x, float y, float z, float sx, float sy, float sz, float point_size, float angle)
{
  float dist, psize, scale;
  scale = sx;
  dist = sqrt((x - playerx) * (x - playerx) + (y - playery) * (y - playery) + (z - playerz) * (z - playerz));
  psize = SCREEN_W * (scale) / dist;
  if (psize < .01)
    return;
  glPointSize(psize * point_size);
  // printf("Points: %f, dist: %f\n",psize,dist);
  glPushMatrix();
  glScalef(sx, sy, sz);
  glTranslatef(x / sx, y / sy, z / sz);
  glRotatef(angle, 0, 1, 0);
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

void
draw_buffer(BUFFER* b, float x, float y, float z, float scale, float angle)
{
  draw_buffer_ex(b, x, y, z, scale, scale, scale, 1.0f, angle);
}

void
draw_eyes(FURBALL* f)
{
  float x1 = 1.2, y = 10, z1 = 1, x2 = 1.2, z2 = 1, s;
  x1 *= cos(f->a) * (f->ultimate ? 1 : 10);
  z1 *= sin(f->a) * (f->ultimate ? 1 : 8);

  x2 *= cos(-f->a) * (f->ultimate ? 1 : 10);
  z2 *= sin(-f->a) * (f->ultimate ? 1 : 8);
  y = (f->ultimate ? 1.5 : 10);
  glGetFloatv(GL_POINT_SIZE, &s);
  if (f->ultimate)
    {
      glGetFloatv(GL_LINE_WIDTH, &s);
      s *= .25;
    }
  if (s < 2.0f)
    return;
  glPointSize(s * 4);
  glPushMatrix();
  glScalef(f->scale, f->scale, f->scale);
  glTranslatef(f->x / f->scale, f->y / f->scale, f->z / f->scale);

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

void
draw_furball_ultimate(FURBALL* f)
{
  int c, x, y, z, v;
  float *vertex, *basevertex, hair, dir[3], dist, lsize, deviation, *colour, *basecolour, colour_deviation;
  deviation = 1;
  colour_deviation = .2;
  vertex = f->mine->vtx;
  basevertex = f->buf->vtx;
  colour = f->mine->clr;
  basecolour = f->buf->clr;
  for (x = 0; x < f->density; x++)
    {
      for (y = 0; y < f->density; y++)
        {
          for (z = 0; z < f->density; z++)
            {
              for (c = 0; c < TRAIL; c++)
                {
                  // if (!c)
                  //{
                  dir[0] = (f->trail[c][0] - f->x); /// f->scale;
                  dir[1] = (f->trail[c][1] - f->y); /// f->scale;
                  dir[2] = (f->trail[c][2] - f->z); /// f->scale;
                  //} else {
                  // for (v=0;v<3;v++)
                  // {
                  //   dir[v]=f->trail[c][v]-x;
                  // }
                  // }
                  // hair=(c*1.0)/(cuts-1);
                  for (v = 0; v < 3; v++)
                    vertex[v] = basevertex[v] + (dir[v]); //+((((rand()%512)/256.0f)-1.0f)*deviation);
                  // for (v=0;v<3;v++)
                  // {
                  // dir[v]=x-f->trail[c][v];
                  // }
                  dir[0] = f->trail[c + 1][0] - f->x;
                  dir[1] = f->trail[c + 1][1] - f->y;
                  dir[2] = f->trail[c + 1][2] - f->z;
                  for (v = 0; v < 3; v++)
                    vertex[3 + v] = basevertex[3 + v] + (dir[v]); //+((((rand()%512)/256.0f)-1.0f)*deviation);
                  // for (v=0;v<2;v++) colour[v]=basecolour[v]+((((rand()%512)/256.0f)-1.0f)*colour_deviation);
                  vertex += 6;
                  basevertex += 6;
                  colour += 6;
                  basecolour += 6;
                }
            }
        }
    }

  // printf("Drawing furball at %f %f %f\n",dir[0],dir[1],dir[2]);
  dist = sqrt((f->x - playerx) * (f->x - playerx) + (f->y - playery) * (f->y - playery) + (f->z - playerz) * (f->z - playerz));
  lsize = 640 / dist;
  // printf("Drawing furball at %f %f %f : %f - %f\n",f->x,f->y,f->z,lsize,dist);
  if (lsize < .01)
    return;
  glLineWidth(lsize);
  draw_buffer(f->mine, f->x, f->y, f->z, f->scale, 0);
  /*
    for ( c=0;c<TRAIL;c++)
        {
          draw_buffer(f->buf,f->trail[c][0],f->trail[c][1],f->trail[c][2],3);
        }
        draw_buffer(f->buf,f->x,f->y,f->z,3);
        */
}

void
draw_furball_normal(FURBALL* f)
{
  float stretch, displace;

  stretch = 1.0 + pow(ABS(f->y - f->trail[1][1]) * .1, 2);
  displace = f->y - f->trail[1][1];
  // printf("Stretch: %f displace: %f at %08x\n",stretch,displace,f->mine);
  if (stretch > 3.0f)
    stretch = 3.0f;
  // if (displace>3.0f) displace=3.0f;
  // printf
  draw_buffer_ex(f->mine, f->x, f->y - displace, f->z, f->scale, f->scale * stretch, f->scale, 1.0f, 0.0f);
}

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
          // printf("Dying: %i\n",f->dying);
          draw_blood(f->blood, PARTICLES);
        }
    }
}

void
draw_ballz()
{
  int c;
  // printf("drawing ballz!\n");
  for (c = 0; c < BALLZ; c++)
    {
      //  printf("drawing ball %i\n",c);
      draw_furball(ballz[c]);
    }
  glLineWidth(LINE_WIDTH);
}

void
draw_ents()
{
  int c;
  float x, y, z, dist;
  for (c = 0; c < num_ents; c++)
    {
      x = ents[c].x - playerx;
      // y=ents[c].y-playery;
      z = (ents[c].z) - playerz;
      dist = (x * x + z * z) - 32 * 32 * ents[c].s * ents[c].s;
      // printf("Distance of entity %i is %f\n",c,sqrt(dist));
      // printf("Entity %f %f player %f %f\n",ents[c].x,ents[c].z,playerx,playerz);
      if (dist < DRAW_DIST * DRAW_DIST)
        {
          // printf("Drawing entity %i\n",c);
          draw_buffer(ents[c].buf, ents[c].x, ents[c].y, ents[c].z, ents[c].s, ents[c].a);
        }
    }
}

void
draw_tree()
{
  clock_t clk;
  int c, n, x, y, pix, s_x, s_y, e_x, e_y;
  float fx, fy, fs, r, g, b, xx, yy;
  clk = clock();
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
              // else glColor3f(.4,.4,.4);
              glColor3f(r, g, b);

              glVertex3f(fx, 0.0, fy);
              glVertex3f(fx, 0.0, fy + fs);
              glVertex3f(fx + fs, 0.0, fy + fs);
              glVertex3f(fx + fs, 0.0, fy);
            }
        }
    }

  glEnd();
  // printf("Tiles took %u msec\n",clock()-clk);clk = clock();
  glMatrixMode(GL_MODELVIEW);

  draw_ents();
  // printf("Entities took %u msec\n",clock()-clk);clk = clock();
  // draw_buffer(tall_tree,100,0,100,2);

  draw_buffer(grass, 0, 0, 0, 1, 0);
  // printf("Grass took %u msec\n",clock()-clk);clk = clock();
  // draw_furball(fur);
  draw_ballz();
  // printf("Balls took %u msec\n",clock()-clk);clk = clock();
  draw_blood(blood, PARTICLES);
  // printf("Blood took %u msec\n",clock()-clk);clk = clock();
}

void
update_blood(BLOOD* b, int num)
{
  int c, pic_x, pic_y, pr, pg, pb, pix;
  float ax, ay, az;
  for (c = 0; c < num; c++)
    {
      if (b[c].alive)
        {
          b[c].x += b[c].vx;
          b[c].y += b[c].vy;
          b[c].z += b[c].vz;
          b[c].vy -= .6;
          if (b[c].y < 0)
            {
              b[c].alive = 0;
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

BUFFER*
generate_cloud_double(const char* filename_x, const char* filename_z, int density)
{
  BUFFER* b;
  BITMAP *bmpx, *bmpz;
  int x, y, z, d, size_x, size_y, size_z, point_counter, img_x, img_y;
  int img_z, col1, col2, c_x, c_y, c_z;
  float deviation, *vertex, *colour, colour_deviation;
  // printf("Creating cloud from %s...\n",filename);
  deviation = .4f;
  colour_deviation = .03f;
  bmpx = load_bmp(filename_x, 0);
  bmpz = load_bmp(filename_z, 0);
  size_x = bmpx->w;
  size_y = bmpx->h;
  size_z = bmpz->w;

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

                  img_x = x; // floor(sqrt((x-c_x)*(x-c_x)+(z-c_z)*(z-c_z)))+c_x;
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

                  img_x = x; // floor(sqrt((x-c_x)*(x-c_x)+(z-c_z)*(z-c_z)))+c_x;
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

FURBALL*
spawn_furball(float x, float y, float z, BUFFER* b, int density, int ultimate, float red, float green, float blue)
{
  FURBALL* f;
  int bufsize, c;
  float *colour, *basecolour;
  bufsize = density * density * density * TRAIL * (ultimate ? 12 : 3) * sizeof(float);
  f = malloc(sizeof(FURBALL));
  f->buf = b;
  f->mine = malloc(sizeof(BUFFER));
  // printf("Mine is %08x\n",f->mine);
  f->mine->vtx = malloc(bufsize);
  f->mine->clr = malloc(bufsize);
  f->mine->size = b->size;
  f->mine->hardware = 0;
  f->mine->mode = b->mode;
  memcpy(f->mine->vtx, f->buf->vtx, bufsize);
  memcpy(f->mine->clr, f->buf->clr, bufsize);
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

BUFFER*
generate_furball_ultimate(float size, int cuts, int density)
{
  BUFFER* b;
  int x, y, z, c, point_counter;
  float fx, fy, fz, *vertex, *colour;
  float deviation, colour_deviation;
  deviation = size * .1f;
  // deviation=density*.333;
  colour_deviation = .03;
  float hair_colour, hair_distance;
  // printf("Dimensions are: %i x %i x %i\n", size_x,size_y,size_z);
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
              fx = (x * M_PI * 2) / (density - 1);
              fy = (y * M_PI) / (density - 1);
              fz = (z * M_PI * 2) / (density - 1);
              for (c = 0; c < cuts; c++)
                {

                  hair_colour = (.5 + .5 * (c * 1) / (cuts));
                  hair_distance = size + (c * size) / (cuts);
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
                  hair_colour = ((c + 1) * 1.0f) / (cuts);
                  hair_distance = size + ((c + 1) * size) / (cuts);
                  vertex[3] = hair_distance * cos(fx) * sin(fy) + ((((rand() % 512) / 256.0f) - 1.0f) * deviation);
                  vertex[4] = hair_distance * cos(fy) + ((((rand() % 512) / 256.0f) - 1.0f) * deviation);
                  vertex[5] = hair_distance * sin(fz) * sin(fy) + ((((rand() % 512) / 256.0f) - 1.0f) * deviation);
                  colour[3] = hair_colour + ((((rand() % 512) / 256.0f) - 1.0f) * colour_deviation);
                  colour[4] = hair_colour + ((((rand() % 512) / 256.0f) - 1.0f) * colour_deviation);
                  colour[5] = hair_colour + ((((rand() % 512) / 256.0f) - 1.0f) * colour_deviation);
                  vertex += 6;
                  colour += 6;
                  point_counter += 2;
                }
            }
        }
    }
  b->mode = GL_LINES;

  b->size = point_counter;
  // vboize(b);
  printf("Created furball of %i points!\n", point_counter);

  return b;
}

BUFFER*
generate_furball_normal(float size, int cuts, int density)
{
  BUFFER* b;
  int x, y, z, c, point_counter;
  float fx, fy, fz, *vertex, *colour;
  float deviation, colour_deviation, colour_factor;
  deviation = size * .3f;
  // deviation=density*.333;
  colour_deviation = .08;
  float hair_colour, hair_distance;
  // printf("Dimensions are: %i x %i x %i\n", size_x,size_y,size_z);
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
              fx = (x * M_PI * 2) / (density - 1);
              fy = (y * M_PI) / (density - 1);
              fz = (z * M_PI * 2) / (density - 1);
              for (c = 0; c < cuts; c++)
                {

                  vertex[0] = size * cos(fx) * sin(fy) + ((((rand() % 512) / 256.0f) - 1.0f) * deviation);
                  vertex[1] = size * cos(fy) + ((((rand() % 512) / 256.0f) - 1.0f) * deviation);
                  vertex[2] = size * sin(fz) * sin(fy) + ((((rand() % 512) / 256.0f) - 1.0f) * deviation);
                  colour_factor = sqrt(vertex[0] * vertex[0] + vertex[1] * vertex[1] + vertex[2] * vertex[2]) / size;
                  // colour_factor*=colour_factor;
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
  // vboize(b);
  printf("Created furball of %i points!\n", point_counter);

  return b;
}

void
create_ballz()
{
  int c, is_good;
  float x, y, z, r, g, b;
  for (c = 0; c < BALLZ; c++)
    {
      is_good = !(rand() % 8);
      if (LOFI)
        is_good = 0;
      r = ((rand() % 256) / 256.0f);
      g = ((rand() % 256) / 256.0f);
      b = ((rand() % 256) / 256.0f);
      x = ((rand() % 1024) / 1024.0f) * WORLD_SIZE;
      y = HEIGHT;
      z = ((rand() % 1024) / 1024.0f) * WORLD_SIZE;
      ballz[c] = spawn_furball(x, y, z, is_good ? ultimate_furball : casual_furball, 8, is_good, r, g, b);
      ballz[c]->iq = 0;
      ballz[c]->smart = .05f + ((rand() % 256) / 256.0f) * .15;
      if (is_good)
        ballz[c]->smart *= .2;
      if (!is_good)
        ballz[c]->scale = .2f + ((rand() % 256) / 256.0f) * .15;
      else
        ballz[c]->scale = 1.0f + ((rand() % 256) / 256.0f) * 3.0f;
      ballz[c]->dying = 0;
      if (c == 8)
        printf("mine is %08x!\n", ballz[c]->mine);
      // ballz[c]->a= b=((rand()%256)/256.0f);
    }
}

void
update_ballz()
{
  int c, i;
  float x, y, z, r, g, b;
  FURBALL* fur;
  // printf("updateing ballz!\n");
  for (c = 0; c < BALLZ; c++)
    {
      if (ballz[c]->exists)
        {
          //  printf("updateing ball %i !\n",c);
          fur = ballz[c];
          //    if (c==7) printf("mine at %i is %08x!\n",c,ballz[8]->mine);
          fur->iq -= fur->smart;

          if (fur->iq <= 0)
            {
              fur->iq = (IQ + (rand() % IQ)) * M_PI;

              fur->speed = ((rand() % 256) / 256.0f) * SPEED;
              fur->a = ((rand() % 256) / 256.0f) * M_PI * 2;
              fur->bounce_rate = 10.0f + ((rand() % 256) / 256.0f) * 40.0f;
            }
          //      if (c==7) printf("mine at %i is %08x!\n",c,ballz[8]->mine);
          fur->bounce += .1f * SPEED;
          fur->x += cos(fur->a) * fur->speed;
          fur->z += sin(fur->a) * fur->speed;
          fur->y = fur->scale * (fur->ultimate ? 1 : 8) + ABS(sin(fur->iq) * fur->bounce_rate);
          if (fur->x < 0)
            fur->x = 0;
          if (fur->z < 0)
            fur->z = 0;
          if (fur->x > WORLD_SIZE)
            fur->x = WORLD_SIZE;
          if (fur->z > WORLD_SIZE)
            fur->z = WORLD_SIZE;
          // if (c==7) printf("mine at %i is %08x!\n",c,ballz[8]->mine);
          for (i = TRAIL; i > 0; i--)
            {

              fur->trail[i][0] = fur->trail[i - 1][0];
              fur->trail[i][1] = fur->trail[i - 1][1];
              fur->trail[i][2] = fur->trail[i - 1][2];
            }
          // if (c==7) printf("mine at %i is %08x!\n",c,ballz[8]->mine);
          // fur->y=playery;
          // fur->x=playerz;
          fur->trail[0][0] = fur->x;
          fur->trail[0][1] = fur->y;
          fur->trail[0][2] = fur->z;
        }
      else if (ballz[c]->dying)
        {
          // printf("%i is dying %i\n",c,ballz[c]->dying);
          update_blood(ballz[c]->blood, PARTICLES);
          ballz[c]->dying--;
        }
      // if (c==7) printf("mine at %i is %08x!\n",c,ballz[8]->mine);haha
    }
}

BUFFER*
generate_world_map(const char* grass_file, const char* height_file, const char* entity_file, int density)
{
  BUFFER* b;
  BITMAP *bmpg, *bmph, *bmph_t, *bmpg_t;
  int x, y, z, d, size_x, size_y, line_counter, img_x, img_y, col;
  float deviation, *vertex, *colour, colour_deviation, grass_x, grass_y, grass_h, patch_size;
  printf("Creating grass from %s...\n", grass_file);
  deviation = .4f;
  colour_deviation = .1f;
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

  b = malloc(sizeof(BUFFER));
  b->vtx = malloc(size_x * size_y * density * 6 * sizeof(float));
  b->clr = malloc(size_x * size_y * density * 6 * sizeof(float));
  b->hardware = 0;
  line_counter = 0;
  vertex = b->vtx;
  colour = b->clr;
  patch_size = WORLD_SIZE / bmpg->w;
  printf("Beggining clouding...\n");
  for (x = 0; x < size_x; x++)
    {
      for (y = 0; y < size_y; y++)
        {
          for (d = 0; d < density; d++)
            {
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
  bmph = load_bmp(entity_file, 0);
  size_x = bmph->w;
  size_y = bmph->h;
  for (x = 0; x < size_x; x++)
    {
      for (y = 0; y < size_y; y++)
        {
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
  // fur = spawn_furball(20,0,20,casual_furball,8,0,.4,.5,.6);
  grass = generate_world_map("grassmap.bmp", "grassheight.bmp", "entitymap.bmp", LOFI ? 2 : 4);
  create_ballz();
  // world[0] = generate_cloud_single("horsex.bmp",2);
}

void
update_furballs()
{
  int c;
}

void
timer_proc(void)
{
  int c;
  float move_spd = MOVE_SPEED;
  clock_t clk;
  clk = clock();
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
      if (talking)
        {
          talk_counter--;
          if (talk_counter <= 0)
            {
              talk_counter = TALK_DELAY + (rand() % TALK_DELAY);
              // play_sample(furtalk[(rand()%16)>>1],128,128,1000,0);
            }
        }
      // printf("mine inb4 timer is %08x!\n",ballz[8]->mine);
      if (whistle_timeout)
        whistle_timeout--;
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
      if (key[KEY_LSHIFT])
        move_spd *= 2;
      get_mouse_mickeys(&mx, &my);
      position_mouse(SCREEN_W / 2, mouse_y);
      update_blood(blood, PARTICLES);
      lookx = (M_PI / 2) * (1.0f * (mouse_y - (SCREEN_W / 2))) / (SCREEN_W / 2.0f);
      looky += (M_PI) * (mx * 1.0f) / (SCREEN_W / 2.0f);
      // if (key[KEY_L])
      //{
      // printf("pos: %f %f %f look: %f %f %f",playerx,playery,playerz,lookx,looky,lookz);
      //}
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
      // if (key[KEY_Q]) explode(blood,PARTICLES,10,10,10);
      if (!key[KEY_UP] && !key[KEY_DOWN] && state != IDLE)
        state = IDLE;
      if ((playery - HEIGHT) < .2)
        {
          // printf("b");
          // if (key[KEY_LCONTROL])
          // {
          //     bounce*=1.7;
          // }
          //  else bounce*=.3;
        }
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
      if (mouse_b)
        shoot();
    }
  bounce = BOUNCE_RATE;
  // if (bounce>400) bounce=400;
  // printf("%f\n",bounce);
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
  /*
      if (key[KEY_SPACE])
      {
          if (lookx<M_PI/4) lookx+=.1;
      }
      else
      {
          if (lookx>0.0) lookx-=.2;
      }
      */
  if (playerx > (WORLD_SIZE - 50.0))
    playerx = (WORLD_SIZE - 50.0);
  if (playerz > (WORLD_SIZE - 50.0))
    playerz = (WORLD_SIZE - 50.0);
  if (playerx < 50.0)
    playerx = 50.0;
  if (playerz < 50.0)
    playerz = 50.0;

  // if (key[KEY_A]) gameover=1;
  //  printf("mine inb4 update is %08x!\n",ballz[8]->mine);

  update_ballz();
  // printf("Update took %u msec\n",clock()-clk);
}

void
draw(void)
{
  clock_t clk;
  int c;
  float r, g, b, mx, my, bh;
  float bg, bg2, ns;
  char nums[16];
  int left = 0;
  for (c = 0; c < BALLZ; c++)
    {
      if (ballz[c]->exists)
        left++;
    }
  if (key[KEY_Y])
    left = 0;
  if (!left)
    game_state = OUTRO;
  if (left < 100)
    talking = 0;
  sprintf(nums, "%03i:%03i", left, BALLZ);
  for (c = 0; c < 7; c++)
    {
      nums[c] -= 0x30;
      //  printf("%i:%i ",c,nums[c]);
    }
  // printf("\n");
  clk = clock();
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glPushMatrix();
  glFrustum(-(4.0 / 3.0), 4.0 / 3.0, -1.0, 1.0, 1.0, DRAW_DIST);
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
  // glRotatef (DEG(lookx), 1, 0, 0);

  glLoadIdentity();
  glTranslatef(0, sin(lookx * 1.333), 0);
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glDisable(GL_DEPTH_TEST);
  glDepthMask(GL_FALSE);
  glBegin(GL_QUADS);

  glColor3f(1.0f, .8, .4);
  glVertex2f(-1, .5);
  glVertex2f(1, .5);
  glColor3f(1.0f, .85, .6);
  glVertex2f(1, .2);
  glVertex2f(-1, .2);

  // glColor3f(1.0f,.8,.4);
  glVertex2f(-1, .2);
  glVertex2f(1, .2);
  glColor3f(1.0f, .9, .8);
  glVertex2f(1, 0);
  glVertex2f(-1, 0);

  // glColor3f(1.0f,1,1);
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
  draw_tree();

  // draw_trees();
  // glFlush();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  // glPointSize(3.0);
  glDisable(GL_TEXTURE_2D);
  glColor3f(1.0, 1.0, 1.0);
  mx = (mouse_x - 160) / 160.0;
  my = -(mouse_y - 120) / 120.0;

  glPointSize(3.0);
  glDisable(GL_DEPTH_TEST);
  if (game_state == PLAY)
    {
      glBegin(GL_POINTS);
      glColor4f(1.0, 1.0, 1.0, 1.0);
      glVertex2f(0, 0);

      glEnd();
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
  allegro_gl_flip();
  // printf("Draw took %u msec\n",clock()-clk);
}

void
timer(void)
{
  tim++;
}
END_OF_FUNCTION(timer);

int
main(int argc, char** argv)
{
  int c;
  float fogc[] = {1.0f, .8, .4, .0};
  char samp[64];
  allegro_init();
  // printf("WAT?\n");
  install_allegro_gl();

  allegro_gl_clear_settings();
  allegro_gl_set(AGL_COLOR_DEPTH, 32);
  allegro_gl_set(AGL_Z_DEPTH, 24);
  allegro_gl_set(AGL_WINDOWED, TRUE);
  allegro_gl_set(AGL_DOUBLEBUFFER, 1);
  allegro_gl_set(AGL_SUGGEST, AGL_COLOR_DEPTH | AGL_Z_DEPTH
                                | AGL_DOUBLEBUFFER | AGL_WINDOWED);

  if (set_gfx_mode(GFX_OPENGL, 640, 480, 0, 0) < 0)
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
  //	show_mouse(screen);

  LOCK_FUNCTION(timer);
  LOCK_VARIABLE(tim);

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
  install_int_ex(timer, BPS_TO_TIMER(30));
  srand(time(NULL));
  // printf("WAT?\n");
  glEnable(GL_TEXTURE_2D);
  generate_stuff();
  // printf("mine in main is %08x!\n",ballz[8]->mine);
  glClearColor(1.0f, .8, .4, .0);

  glEnable(GL_BLEND);
  glEnable(GL_ALPHA_TEST);
  glAlphaFunc(GL_GREATER, .3);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  playerx = 20.0;
  playerz = 20.0;
  looky = M_PI / 4;
  // printf("mine inb4 loop is %08x!\n",ballz[8]->mine);
  get_mouse_mickeys(&mx, &my);

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
      sprintf(samp, "fur%i.wav\0", c + 1);
      printf("Loading %s\n", samp);
      furtalk[c] = load_wav(samp);
      sprintf(samp, "die%i.wav\0", c + 1);
      printf("Loading %s\n", samp);
      die[c] = load_wav(samp);
    }
  playerx = 150;
  playery = HEIGHT;
  playerz = 150;
  lookx = (M_PI / 2) * (1.0f * ((SCREEN_H / 2) - (SCREEN_W / 2))) / (SCREEN_W / 2.0f);
  ;
  looky = 1.639519;
  lookz = 0;
  play_sample(music, 255, 128, 1000, 1);
  do
    {
      // if (tim<=0)
      while (tim > 0)
        {
          timer_proc();
          tim--;
        }
      draw();
      // Sleep(5);
    }
  while (!GetAsyncKeyState(VK_ESCAPE));

  set_gfx_mode(GFX_TEXT, 0, 0, 0, 0);

  return 0;
}
END_OF_MAIN()
