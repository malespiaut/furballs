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
#include <float.h> // for FLT_EPSILON
#include <math.h>  // just math
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>  // standard
#include <string.h> // for isExtensionSupported

#include <SOIL/SOIL.h> // Simple OpenGL Image Loader
#include <alleggl.h>   // gl helper for allegro
#include <allegro.h>   // allegro - a game programming library
////////////////////////////////////////////////////
//////////////////////////////////////////////////// DEFINES
////////////////////////////////////////////////////

#define kScreenWidth 1920
#define kScreenHeight 1200

#define kPi 3.14159265358979323846f               /* pi */
#define kPiDividedBy2 1.57079632679489661923f     /* pi/2 */
#define kPiDividedBy4 0.78539816339744830962f     /* pi/4 */
#define kPiMultipliedBy2 6.28318530717958647692f  /* pi*2 */
#define kPiMultipliedBy4 12.56637061435917295384f /* pi*4 */

// tweakable defines
#define LOFI 0             // this disables VBOs, reduces drawing distance and furball complexity
                           // HiFi version runs smoothly on a modern PC or decent laptop
                           // LoFi version runs smoothly on a 10yo PC
#define kMoveSpeed 1.0f    // camera/player move speed
#define kWorldSize 2000.0f // side of scene square
#define kHaystack 100.0f   // player height above the ground modifier

// defines for calculation simplification
#define FRAND(x) (((xrnd() % (int)x) * 1024.0f) / 1024.0f)
#define DEG(n) ((n)*180.0f / kPi)

// player states
#define WALK 1
#define IDLE 2
#define STOP 3

#define kHeight 10.0f    // height of camera/player
#define kBounceRate 2.0f // player jump rate
#define kNumSize 0.1f    // size of onscreen numbers

// game state
#define INTRO 1
#define PLAY 2
#define OUTRO 3

#define kIQ 2       // furball bounces between change of direction
#define kSpeed 1.5f // furball speed

// graphics quality settings
#if LOFI == 1
#define kBallz 666       // number of furballs
#define kDrawDist 500.0f // rendering / fog distance
#define kLineWidth 24.0f // line width for grass
#define kParticles 128   // number of max onscreen particles
#define kTrail 1         // hair length for hairy furballs
#else
#define kBallz 666
#define kDrawDist 2000.0f
#define kLineWidth 12.0f
#define kParticles 256
#define kTrail 4
#endif

#define kMaxEntities 1024 // max scene entities
#define kSamps 8          // number of furball voice samples
#define kTalkDelay 100    // delay between furball sample playback

////////////////////////////////////////////////////
//////////////////////////////////////////////////// TYPEDEFS
////////////////////////////////////////////////////

typedef struct Vec3 Vec3;
struct Vec3
{
  float x;
  float y;
  float z;
};

typedef struct Vec3i Vec3i;
struct Vec3i
{
  int32_t x;
  int32_t y;
  int32_t z;
};

/*
typedef struct Color Color;
struct Color
{
  float r;
  float g;
  float b;
  float a; // Alpha component for transparency
};
*/

/*
typedef struct Rectangle Rectangle;
struct Rectangle
{
  int32_t left;
  int32_t top;
  int32_t right;
  int32_t bottom;
};
*/

// basic vertex buffer struct
typedef struct Buffer Buffer;
struct Buffer
{
  float* vtx;            // vertex and colour buffers
  float* clr;            // vertex and colour buffers
  uint32_t vtx_handle;   // VBO handles for the above
  uint32_t clr_handle;   // VBO handles for the above
  int32_t size;          // number of vertices
  int32_t hardware;      // is it VBO or array
  int32_t frames;        // number of frames
  int32_t current_frame; // currently playing frame
  uint32_t mode;         // lines or points
};

// blood particle struct
typedef struct Blood Blood;
struct Blood
{
  float x;    // paritcle position
  float y;    // paritcle position
  float z;    // paritcle position
  float r;    // colour
  float g;    // colour
  float b;    // colour
  float a;    // colour
  float vx;   // velocity
  float vy;   // velocity
  float vz;   // velocity
  float s;    // size
  bool alive; // update flag
};

// furball (a dynamic vertex buffer + AI)
typedef struct Furball Furball;
struct Furball
{
  Buffer* buf;                // base vertex buffer for augmentation
  Buffer* mine;               // actually drawn buffer
  float x;                    // position
  float y;                    // position
  float z;                    // position
  float scale;                // size
  float r;                    // colour
  float g;                    // colour
  float b;                    // colour
  bool ultimate;              // is it dynamic (the long-haired one
  float a;                    // angle of movement
  float bounce;               // bounce rate
  float speed;                // speed of movement
  float bounce_rate;          // base bounce rate
  float trail[kTrail + 1][3]; // 'hair' trail (only ultimate)
  size_t density;             // points per voxel
  int32_t exists;             // render flag
  int32_t dying;              // blood render flag
  Blood blood[kParticles];    // blood particles
  float iq;                   // directiona change interval
  float smart;                // iq base value
};

// a buffer with position
typedef struct Entity Entity;
struct Entity
{
  float x;     // position
  float y;     // position
  float z;     // position
  float s;     // size
  float a;     // yaw (angle)
  Buffer* buf; // vertex buffer
};

////////////////////////////////////////////////////
//////////////////////////////////////////////////// GLOBALS
////////////////////////////////////////////////////

int32_t frames = 0;   // frames rendered
int32_t state = IDLE; // player state
volatile int32_t tim; // threaded timing variable

Vec3 player = {.x = 0.0f, .y = 0.0f, .z = 0.0f};
float lookx = 0.0f; // player look direction
float looky = 0.0f; // player look direction
float lookz = 0.0f; // player look direction
float lookf = 0.0f; // player look direction

Buffer* ultimate_furball = NULL; // base buffer for ultimate furball (with ahir)
Buffer* casual_furball = NULL;   // base buffer for simple furball

// vertex buffers for world entities
Buffer* g_grass = NULL;
Buffer* tall_tree = NULL;
Buffer* withered_bush = NULL;
Buffer* cactus = NULL;
Buffer* palm = NULL;
Buffer* stick = NULL;
Buffer* pine = NULL;
Buffer* hut = NULL;
Buffer* fence = NULL;
Buffer* church = NULL;
Buffer* brickhouse = NULL;

BITMAP* ground_bmp = NULL;       // ground colour
float bounce = kBounceRate;      // player bounce rate
Furball* ballz[kBallz] = {NULL}; // furballs
int32_t mx = 0;                  // mouse delta position
int32_t my = 0;                  // mouse delta position
Blood g_blood[kParticles] = {0}; // particles
Entity ents[kMaxEntities] = {0}; // entities (map objects) array
size_t num_ents = 0;             // number of entities
int32_t game_state = INTRO;      // game state
uint32_t numbers = 0;            // textures
uint32_t intro = 0;              // textures
uint32_t outro = 0;              // textures

SAMPLE* youknow = NULL;
SAMPLE* furtalk[kSamps] = {NULL};
SAMPLE* die[kSamps] = {NULL};
SAMPLE* slash = NULL;
SAMPLE* music = NULL;
SAMPLE* shot = NULL;
SAMPLE* whistle = NULL; // sounds

int32_t talk_counter = kTalkDelay;
int32_t talking = 1;           // furball talk timers
int32_t whistle_timeout = 100; // furball whistle timeout

bool game_over = false;

////////////////////////////////////////////////////
//////////////////////////////////////////////////// HELPER FUNCTIONS
////////////////////////////////////////////////////

// checks if GL extension is supported (stolen)
static int
isExtensionSupported(const char* extension)
{
  /* Extension names should not have spaces. */
  uint8_t* where = (uint8_t*)strchr(extension, ' ');
  if (where || *extension == '\0')
  {
    return 0;
  }
  const uint8_t* extensions = glGetString(GL_EXTENSIONS);
  /* It takes a bit of care to be fool-proof about parsing the
     OpenGL extensions string. Don't be fooled by sub-strings,
     etc. */
  const uint8_t* start = extensions;
  for (;;)
  {
    where = (uint8_t*)strstr((const char*)start, extension);
    if (!where)
    {
      break;
    }
    uint8_t* terminator = where + strlen(extension);
    if (where == start || *(where - 1) == ' ')
    {
      if (*terminator == ' ' || *terminator == '\0')
      {
        return 1;
      }
    }
    start = terminator;
  }
  return 0;
}

// xoshiro256+

// xoshiro256+ (a single plus sign) is used here,
// rather than xoshiro256++ (two plus signs), because
// it is described by its designers as the “best and
// fastest generator for floating-point numbers”

uint64_t rngstate[4] = {0};

static inline uint64_t
rotl(const uint64_t x, int k)
{
  return (x << k) | (x >> (64 - k));
}

// using an implementation of splitmix64, as recommended
// by the authors of xoshiro256+, to fill the rngstate array.
static void
xrnd_seed(uint64_t seed)
{
  rngstate[0] = seed;
  for (size_t i = 1; i < 4; ++i)
  {
    rngstate[i] = (rngstate[i - 1] += 0x9e3779b97f4a7c15);
    rngstate[i] = (rngstate[i] ^ (rngstate[i] >> 30)) * 0xbf58476d1ce4e5b9;
    rngstate[i] = (rngstate[i] ^ (rngstate[i] >> 27)) * 0x94d049bb133111eb;
    rngstate[i] ^ (rngstate[i] >> 31);
  }
}

static uint64_t
xrnd(void)
{
  const uint64_t result = rngstate[0] + rngstate[3];

  const uint64_t t = rngstate[1] << 17;

  rngstate[2] ^= rngstate[0];
  rngstate[3] ^= rngstate[1];
  rngstate[1] ^= rngstate[2];
  rngstate[0] ^= rngstate[3];

  rngstate[2] ^= t;

  rngstate[3] = rotl(rngstate[3], 45);

  return result;
}

////////////////////////////////////////////////////
//////////////////////////////////////////////////// MATH FUNCTIONS
////////////////////////////////////////////////////

// vector length
static float
length3v(const float* a)
{
  return sqrtf(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

// 2 point distance
static float
vec3_dist(Vec3 a, Vec3 b)
{
  float dx = b.x - a.x;
  float dy = b.y - a.y;
  float dz = b.z - a.z;
  return sqrtf(dx * dx + dy * dy + dz + dz);
}

// dot product
static float
vec3_dot(Vec3 a, Vec3 b)
{
  return (a.x * b.x + a.y * b.y + a.z * b.z);
}

////////////////////////////////////////////////////
//////////////////////////////////////////////////// RENDERING FUNCTIONS
////////////////////////////////////////////////////

// creates VBO from vertex array
static void
vboize(Buffer* b)
{
  if (LOFI)
  {
    return;
  }
  if (!isExtensionSupported("GL_ARB_vertex_buffer_object"))
  {
    return;
  }

  ptrdiff_t nb_elements = 0;
  if (b->size > 0)
  {
    nb_elements = (ptrdiff_t)sizeof(float) * 3 * b->size;
  }

  b->hardware = 1;
  glGenBuffers(1, &(b->vtx_handle));
  glBindBuffer(GL_ARRAY_BUFFER, b->vtx_handle);
  glBufferData(GL_ARRAY_BUFFER, nb_elements, b->vtx, GL_STATIC_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  glGenBuffers(1, &(b->clr_handle));
  glBindBuffer(GL_ARRAY_BUFFER, b->clr_handle);
  glBufferData(GL_ARRAY_BUFFER, nb_elements, b->clr, GL_STATIC_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  free(b->vtx);
  free(b->clr);
}

// renders a Buffer at desired position
static void
draw_buffer_ex(Buffer* b, float x, float y, float z, float sx, float sy, float sz, float point_size, float angle)
{
  float dist = sqrtf((x - player.x) * (x - player.x) + (y - player.y) * (y - player.y) + (z - player.z) * (z - player.z));

  // calculates point size from distance
  // sx is a scaling factor
  float psize = (float)kScreenWidth * sx / dist;
  if (psize - 0.01f < FLT_EPSILON)
  {
    return;
  }

  glPointSize(psize * point_size);
  glPushMatrix();
  glScalef(sx, sy, sz);
  glTranslatef(x / sx, y / sy, z / sz);
  glRotatef(angle, 0.0f, 1.0f, 0.0f);

  // draws a VBO or vertex array
  if (b->hardware)
  {
    glBindBuffer(GL_ARRAY_BUFFER, b->vtx_handle);
    glVertexPointer(3, GL_FLOAT, 0, NULL);
    glBindBuffer(GL_ARRAY_BUFFER, b->clr_handle);
    glColorPointer(3, GL_FLOAT, 0, NULL);
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
static void
draw_buffer(Buffer* b, float x, float y, float z, float scale, float angle)
{
  draw_buffer_ex(b, x, y, z, scale, scale, scale, 1.0f, angle);
}

// draws eyes for a furball
static void
draw_eyes(Furball* f)
{
  float x1 = 1.2f;
  float y = 10.0f;
  float z1 = 1.0f;
  float x2 = 1.2f;
  float z2 = 1.0f;
  float s = 0.0f;

  // haired furballs are bigger
  x1 *= cosf(f->a) * (f->ultimate ? 1.0f : 10.0f);
  z1 *= sinf(f->a) * (f->ultimate ? 1.0f : 8.0f);

  x2 *= cosf(-f->a) * (f->ultimate ? 1.0f : 10.0f);
  z2 *= sinf(-f->a) * (f->ultimate ? 1.0f : 8.0f);
  y = (f->ultimate ? 1.5f : 10.0f);

  // gets point size for this furball
  glGetFloatv(GL_POINT_SIZE, &s);
  if (f->ultimate)
  {
    // gets line width if it's a hairy one
    glGetFloatv(GL_LINE_WIDTH, &s);
    s *= 0.25f;
  }
  // if the size is too small, no use drawing eyes
  if (s < 2.0f)
  {
    return;
  }
  glPointSize(s * 4.0f);
  glPushMatrix();
  glScalef(f->scale, f->scale, f->scale);
  glTranslatef(f->x / f->scale, f->y / f->scale, f->z / f->scale);

  // draws eyes
  glBegin(GL_POINTS);
  glColor3f(1.0f, 1.0f, 1.0f);
  glVertex3f(x1, y, z1);
  glVertex3f(x2, y, z2);
  glEnd();
  glPointSize(s * 2.0f);
  glBegin(GL_POINTS);
  glColor3f(0.0f, 0.0f, 0.0f);
  glVertex3f(x1 * 1.2f, y, z1 * 1.2f);
  glVertex3f(x2 * 1.2f, y, z2 * 1.2f);
  glEnd();
  glPopMatrix();
}

// draws a long-haired furball
// pretty much more CPU intensive than it seems
// thus ultimate furballs are disabled in LOFI mode
static void
draw_furball_ultimate(Furball* furball)
{
  float dir[3] = {0.0f};

  // fetches base vertex buffer for augmentation
  float* vertex = furball->mine->vtx;
  float* basevertex = furball->buf->vtx;

  // augments the buffer
  for (size_t x = 0; x < furball->density; ++x)
  {
    for (size_t y = 0; y < furball->density; ++y)
    {
      for (size_t z = 0; z < furball->density; ++z)
      {
        for (size_t i = 0; i < kTrail; ++i)
        {
          dir[0] = (furball->trail[i][0] - furball->x);
          dir[1] = (furball->trail[i][1] - furball->y);
          dir[2] = (furball->trail[i][2] - furball->z);
          for (size_t v = 0; v < 3; ++v)
          {
            vertex[v] = basevertex[v] + (dir[v]);
          }

          dir[0] = furball->trail[i + 1][0] - furball->x;
          dir[1] = furball->trail[i + 1][1] - furball->y;
          dir[2] = furball->trail[i + 1][2] - furball->z;
          for (size_t v = 0; v < 3; ++v)
          {
            vertex[3 + v] = basevertex[3 + v] + (dir[v]);
          }
          vertex += 6;
          basevertex += 6;
        }
      }
    }
  }

  // calculates line width and draws

  float dist = sqrtf((furball->x - player.x) * (furball->x - player.x) + (furball->y - player.y) * (furball->y - player.y) + (furball->z - player.z) * (furball->z - player.z));
  float lsize = kScreenWidth / dist;
  // if it's too small, no use drawing
  if (lsize < 0.01f)
  {
    return;
  }
  glLineWidth(lsize);
  draw_buffer(furball->mine, furball->x, furball->y, furball->z, furball->scale, 0);
}

// draws a simple furball from a vertex buffer
// stretches it a bit according to vertical velocity
// and displaces a bit (actually eyes look displaced this way)
static void
draw_furball_normal(Furball* furball)
{
  float stretch = 1.0f + powf(ABS(furball->y - furball->trail[1][1]) * 0.1f, 2.0f);
  float displace = furball->y - furball->trail[1][1];
  if (stretch - 3.0f > FLT_EPSILON)
  {
    stretch = 3.0f;
  }

  draw_buffer_ex(furball->mine, furball->x, furball->y - displace, furball->z, furball->scale, furball->scale * stretch, furball->scale, 1.0f, 0.0f);
}

// draws a blood particles array and makes coffee
static void
draw_blood(Blood* blood, size_t num)
{
  glPointSize(20.0f);
  glBegin(GL_POINTS);
  for (size_t i = 0; i < num; ++i)
  {
    glColor3fv(&(blood[i].r));
    glVertex3fv(&(blood[i].x));
  }
  glEnd();
}

// draws a furball (picks one of the above and checks distance)
static void
draw_furball(Furball* furball)
{
  Vec3 position = {
    .x = furball->x - player.x,
    .y = furball->y - player.y,
    .z = furball->z - player.z};

  // if ((position.x * position.x + position.y * position.y + position.z * position.z) < kDrawDist * kDrawDist)
  if (vec3_dot(position, position) < kDrawDist * kDrawDist)
  {
    if (furball->exists)
    {
      if (furball->ultimate)
      {
        draw_furball_ultimate(furball);
      }
      else
      {
        draw_furball_normal(furball);
      }
      draw_eyes(furball);
    }
    if (furball->dying)
    {
      draw_blood(furball->blood, kParticles);
    }
  }
}

// draws all furballs
static void
draw_ballz(void)
{
  for (size_t i = 0; i < kBallz; ++i)
  {
    draw_furball(ballz[i]);
  }
  glLineWidth(kLineWidth);
}

// draws all entities
// actually draws only the ones that are close enough
static void
draw_ents(void)
{
  for (size_t i = 0; i < num_ents; ++i)
  {
    float x = ents[i].x - player.x;
    float z = (ents[i].z) - player.z;
    float dist = (x * x + z * z) - 32.0f * 32.0f * ents[i].s * ents[i].s;
    if (dist < kDrawDist * kDrawDist)
    {
      draw_buffer(ents[i].buf, ents[i].x, ents[i].y, ents[i].z, ents[i].s, ents[i].a);
    }
  }
}

// apparently, this draws ground.
// ground changes colour dynamically, and is not a textured quad
// but a series of quads, drawn only near the player
// this saves cpu and gpu from dynamic texture updates
// it also draws everything else :P
static void
draw_tree(void)
{
  // draws array of ground quads
  glBegin(GL_QUADS);
  for (int32_t x = 0; x < ground_bmp->w; ++x)
  {
    for (int32_t y = 0; y < ground_bmp->h; ++y)
    {
      float fx = ((float)x * kWorldSize) / (float)ground_bmp->w;
      float fy = ((float)y * kWorldSize) / (float)ground_bmp->w;
      float fs = kWorldSize / (float)ground_bmp->w;
      float xx = fx - player.x;
      float yy = fy - player.z;

      // if close enough, draws
      if (xx * xx + yy * yy < kDrawDist * kDrawDist)
      {
        int32_t c = (x + y) & 1;
        int32_t pix = getpixel(ground_bmp, x, y);
        float r = (float)getr(pix) / 256.0f;
        float g = (float)getg(pix) / 256.0f;
        float b = (float)getb(pix) / 256.0f;
        if (c)
        {
          r *= 0.8f;
          g *= 0.8f;
          b *= 0.8f;
        }

        glColor3f(r, g, b);
        glVertex3f(fx, 0.0f, fy);
        glVertex3f(fx, 0.0f, fy + fs);
        glVertex3f(fx + fs, 0.0f, fy + fs);
        glVertex3f(fx + fs, 0.0f, fy);
      }
    }
  }
  glEnd();
  glMatrixMode(GL_MODELVIEW);

  // draws map object
  draw_ents();

  // draws grass
  draw_buffer(g_grass, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);

  // draws furballs
  draw_ballz();

  // and particles
  draw_blood(g_blood, kParticles);
}

////////////////////////////////////////////////////
//////////////////////////////////////////////////// GENERATOR FUNCTIONS
////////////////////////////////////////////////////

// generates a mesh from one bitmap
// by 'spinning it' around Y axis
static Buffer*
generate_cloud_single(const char* filename, size_t density)
{
  if (!filename)
  {
    return NULL;
  }

  float deviation = 0.4f;
  float colour_deviation = 0.03f;

  BITMAP* bmp = load_bmp(filename, NULL);

  Vec3i size = {.x = bmp->w, .y = bmp->h, .z = bmp->w};
  if (size.x < 0)
  {
    size.x = 0;
  }
  if (size.y < 0)
  {
    size.y = 0;
  }
  if (size.z < 0)
  {
    size.z = 0;
  }

  // allocates a too large buffer, or too small, if of higher density
  Buffer* b = malloc(sizeof(*b));
  b->vtx = malloc((size_t)(size.x * size.y * size.z) * 3 * sizeof(*b->vtx));
  b->clr = malloc((size_t)(size.x * size.y * size.z) * 3 * sizeof(*b->clr));
  b->hardware = 0;

  int32_t point_counter = 0;
  float* vertex = b->vtx;
  float* colour = b->clr;
  Vec3 c = {.x = (float)size.x / 2.0f, .y = 0.0f, .z = (float)size.z / 2.0f};

  for (int32_t x = 0; x < size.x; ++x)
  {
    for (int32_t y = 0; y < size.y; ++y)
    {
      for (int32_t z = 0; z < size.z; ++z)
      {
        for (size_t d = 0; d < density; ++d)
        {
          // tests a voxel against image, where y is y
          // and x is distance from centre
          // if ok, sets colour and displaces a bit
          int32_t img_x = (int32_t)floorf(sqrtf(((float)x - c.x) * ((float)x - c.x) + ((float)z - c.z) * ((float)z - c.z)) + c.x);
          int32_t img_y = bmp->h - y - 1;
          int32_t col = getpixel(bmp, img_x, img_y);
          if (col != makecol(255, 0, 255) && img_x < bmp->w)
          {
            vertex[0] = (float)x + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * deviation) - c.x;
            vertex[1] = (float)y + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * deviation);
            vertex[2] = (float)z + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * deviation) - c.z;
            colour[0] = ((float)getr(col) / 256.0f) + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * colour_deviation);
            colour[1] = ((float)getg(col) / 256.0f) + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * colour_deviation);
            colour[2] = ((float)getb(col) / 256.0f) + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * colour_deviation);
            vertex += 3;
            colour += 3;
            ++point_counter;
          }
        }
      }
    }
  }
  b->size = point_counter;
  b->mode = GL_POINTS;
  vboize(b);
  destroy_bitmap(bmp);
  return b;
}

// generates a mesh from two images
// these are axis projections from front and side
static Buffer*
generate_cloud_double(const char* filename_x, const char* filename_z, size_t density)
{
  if (!filename_x || !filename_z)
  {
    return NULL;
  }

  float deviation = 0.4f;
  float colour_deviation = 0.03f;

  BITMAP* bmpx = load_bmp(filename_x, NULL);
  BITMAP* bmpz = load_bmp(filename_z, NULL);

  Vec3i size = {.x = bmpx->w, .y = bmpx->h, .z = bmpz->w};
  if (size.x < 0)
  {
    size.x = 0;
  }
  if (size.y < 0)
  {
    size.y = 0;
  }
  if (size.z < 0)
  {
    size.z = 0;
  }

  // another bad malloc
  Buffer* b = malloc(sizeof(*b));
  b->vtx = malloc((size_t)(size.x * size.y * size.z) * 3 * sizeof(*b->vtx));
  b->clr = malloc((size_t)(size.x * size.y * size.z) * 3 * sizeof(*b->clr));
  b->hardware = 0;

  int32_t point_counter = 0;
  float* vertex = b->vtx;
  float* colour = b->clr;
  Vec3 c = {.x = (float)size.x / 2.0f, .y = 0.0f, .z = (float)size.z / 2.0f};

  for (int32_t x = 0; x < size.x; ++x)
  {
    for (int32_t y = 0; y < size.y; ++y)
    {
      for (int32_t z = 0; z < size.z; ++z)
      {
        for (size_t d = 0; d < density; ++d)
        {
          // this tests the pixel against two projection images
          // and sets the colour by averaging from those
          int32_t img_x = x;
          int32_t img_z = z;
          int32_t img_y = bmpx->h - y - 1;
          int32_t col1 = getpixel(bmpx, img_x, img_y);
          int32_t col2 = getpixel(bmpz, img_z, img_y);
          if (col1 != makecol(255, 0, 255) && col2 != makecol(255, 0, 255) && img_x < bmpx->w && img_z < bmpz->w)
          {
            vertex[0] = (float)x + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * deviation) - c.x;
            vertex[1] = (float)y + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * deviation);
            vertex[2] = (float)z + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * deviation) - c.z;
            colour[0] = ((float)(getr(col1) + getr(col2)) / 512.0f) + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * colour_deviation);
            colour[1] = ((float)(getg(col1) + getg(col2)) / 512.0f) + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * colour_deviation);
            colour[2] = ((float)(getb(col1) + getb(col2)) / 512.0f) + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * colour_deviation);
            vertex += 3;
            colour += 3;
            ++point_counter;
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
  return b;
}

// generates a mesh from 3 projection images (front, side, top)
static Buffer*
generate_cloud_triple(const char* filename_x, const char* filename_z, const char* filename_y, size_t density)
{
  if (!filename_x || !filename_z || !filename_y)
  {
    return NULL;
  }

  float deviation = 0.4f;
  float colour_deviation = 0.03f;

  BITMAP* bmpx = load_bmp(filename_x, NULL);
  BITMAP* bmpz = load_bmp(filename_z, NULL);
  BITMAP* bmpy = load_bmp(filename_y, NULL);

  Vec3i size = {.x = bmpx->w, .y = bmpx->h, .z = bmpz->w};
  if (size.x < 0)
  {
    size.x = 0;
  }
  if (size.y < 0)
  {
    size.y = 0;
  }
  if (size.z < 0)
  {
    size.z = 0;
  }

  // bad malloc Mk.3
  Buffer* b = malloc(sizeof(*b));
  b->vtx = malloc((size_t)(size.x * size.y * size.z) * 3 * sizeof(*b->vtx));
  b->clr = malloc((size_t)(size.x * size.y * size.z) * 3 * sizeof(*b->clr));
  b->hardware = 0;

  int32_t point_counter = 0;
  float* vertex = b->vtx;
  float* colour = b->clr;
  Vec3 c = {.x = (float)size.x / 2.0f, .y = 0.0f, .z = (float)size.z / 2.0f};

  for (int32_t x = 0; x < size.x; ++x)
  {
    for (int32_t y = 0; y < size.y; ++y)
    {
      for (int32_t z = 0; z < size.z; ++z)
      {
        for (size_t d = 0; d < density; ++d)
        {
          // similar to the above, tests against 3 images and
          // averages the colour
          int32_t img_x = x;
          int32_t img_z = z;
          int32_t img_y = bmpx->h - y - 1;
          int32_t col1 = getpixel(bmpx, img_x, img_y);
          int32_t col2 = getpixel(bmpz, img_z, img_y);
          int32_t col3 = getpixel(bmpy, img_x, img_z);
          if (col3 != makecol(255, 0, 255) && col1 != makecol(255, 0, 255) && col2 != makecol(255, 0, 255) && img_x < bmpx->w && img_z < bmpz->w)
          {
            vertex[0] = (float)x + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * deviation) - c.x;
            vertex[1] = (float)y + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * deviation);
            vertex[2] = (float)z + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * deviation) - c.z;
            colour[0] = ((float)(getr(col1) + getr(col2)) / 512.0f) + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * colour_deviation);
            colour[1] = ((float)(getg(col1) + getg(col2)) / 512.0f) + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * colour_deviation);
            colour[2] = ((float)(getb(col1) + getb(col2)) / 512.0f) + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * colour_deviation);
            vertex += 3;
            colour += 3;
            ++point_counter;
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
  return b;
}

//  this creates a furball instance
static Furball*
spawn_furball(float x, float y, float z, Buffer* b, size_t density, bool ultimate, float red, float green, float blue)
{
  // allocates memory
  Furball* furball = malloc(sizeof(*furball));
  furball->mine = malloc(sizeof(*furball->mine));

  // number of mesh vertices is specified through density
  // it represents, slices, sectors and depth
  size_t bufsize = density * density * density * kTrail * (ultimate ? 12 : 3) * sizeof(*furball->mine->vtx);

  furball->mine->vtx = malloc(bufsize);
  furball->mine->clr = malloc(bufsize);
  furball->buf = b;
  furball->mine->size = b->size;
  furball->mine->mode = b->mode;
  furball->mine->hardware = 0;

  // copies the buffer from base for dynamic updates
  memcpy(furball->mine->vtx, furball->buf->vtx, bufsize);
  memcpy(furball->mine->clr, furball->buf->clr, bufsize);

  // sets colour
  float* colour = furball->mine->clr;
  float* basecolour = furball->buf->clr;
  for (int32_t i = 0; i < furball->mine->size; ++i)
  {
    colour[0] = basecolour[0] * red;
    colour[1] = basecolour[1] * green;
    colour[2] = basecolour[2] * blue;
    colour += 3;
    basecolour += 3;
  }

  // sets position and other properties
  furball->x = x;
  furball->y = y;
  furball->z = z;
  furball->r = red;
  furball->g = green;
  furball->b = blue;
  furball->density = density;
  furball->scale = 3.0f;
  furball->ultimate = ultimate;
  furball->exists = 1;
  return furball;
}

// generates a base buffer for hairy furball
static Buffer*
generate_furball_ultimate(float size, size_t cuts, size_t density)
{
  float deviation = size * 0.1f;
  float colour_deviation = 0.03f;

  // allocates memory
  Buffer* b = malloc(sizeof(*b));
  b->vtx = malloc(density * density * density * cuts * 12 * sizeof(*b->vtx));
  b->clr = malloc(density * density * density * cuts * 12 * sizeof(*b->clr));
  b->hardware = 0;
  int32_t point_counter = 0;
  float* vertex = b->vtx;
  float* colour = b->clr;

  for (size_t x = 0; x < density; ++x)
  {
    for (size_t y = 0; y < density; ++y)
    {
      for (size_t z = 0; z < density; ++z)
      {
        // generates a sphere of vertices
        float fx = ((float)x * kPiMultipliedBy2) / (float)(density - 1);
        float fy = ((float)y * kPi) / (float)(density - 1);
        float fz = ((float)z * kPiMultipliedBy2) / (float)(density - 1);

        // generates protruding hair
        for (size_t i = 0; i < cuts; ++i)
        {

          float hair_colour = 0.5f + 0.5f * (float)i / (float)(cuts);
          float hair_distance = size + ((float)i * size) / (float)(cuts);

          // generates position only for the hair beginning
          // rest is determined by velocity
          // this generates startpoint
          if (i == 0)
          {
            vertex[0] = hair_distance * cosf(fx) * sinf(fy) + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * deviation);
            vertex[1] = hair_distance * cosf(fy) + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * deviation);
            vertex[2] = hair_distance * sinf(fz) * sinf(fy) + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * deviation);
            colour[0] = hair_colour + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * colour_deviation);
            colour[1] = hair_colour + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * colour_deviation);
            colour[2] = hair_colour + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * colour_deviation);
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
          hair_colour = ((float)(i + 1) * 1.0f) / (float)cuts;
          hair_distance = size + ((float)(i + 1) * size) / (float)cuts;
          vertex[3] = hair_distance * cosf(fx) * sinf(fy) + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * deviation);
          vertex[4] = hair_distance * cosf(fy) + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * deviation);
          vertex[5] = hair_distance * sinf(fz) * sinf(fy) + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * deviation);
          colour[3] = hair_colour + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * colour_deviation);
          colour[4] = hair_colour + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * colour_deviation);
          colour[5] = hair_colour + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * colour_deviation);

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

  return b;
}

// egnerates a simple furball
static Buffer*
generate_furball_normal(float size, size_t cuts, size_t density)
{
  float deviation = size * 0.3f;
  float colour_deviation = 0.08f;

  // allocates buffers
  Buffer* b = malloc(sizeof(*b));
  b->vtx = malloc(density * density * density * cuts * 3 * sizeof(*b->vtx));
  b->clr = malloc(density * density * density * cuts * 3 * sizeof(*b->clr));
  b->hardware = 0;

  int32_t point_counter = 0;
  float* vertex = b->vtx;
  float* colour = b->clr;

  for (size_t x = 0; x < density; ++x)
  {
    for (size_t y = 0; y < density; ++y)
    {
      for (size_t z = 0; z < density; ++z)
      {
        // gnerates vertex position in a sphere, and displaces a bit
        float fx = ((float)x * kPiMultipliedBy2) / (float)(density - 1);
        float fy = ((float)y * kPi) / (float)(density - 1);
        float fz = ((float)z * kPiMultipliedBy2) / (float)(density - 1);
        for (size_t i = 0; i < cuts; ++i)
        {
          vertex[0] = size * cosf(fx) * sinf(fy) + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * deviation);
          vertex[1] = size * cosf(fy) + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * deviation);
          vertex[2] = size * sinf(fz) * sinf(fy) + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * deviation);

          float colour_factor = sqrtf(vertex[0] * vertex[0] + vertex[1] * vertex[1] + vertex[2] * vertex[2]) / size;
          colour[0] = colour_factor + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * colour_deviation);
          colour[1] = colour_factor + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * colour_deviation);
          colour[2] = colour_factor + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * colour_deviation);
          vertex += 3;
          colour += 3;
          point_counter++;
        }
      }
    }
  }
  b->mode = GL_POINTS;
  b->size = point_counter;

  return b;
}

// creates all the furballs
static void
create_ballz(void)
{
  for (size_t i = 0; i < kBallz; ++i)
  {
    // if is_good, then it's 'ultimate' == hairy
    bool is_good = !(xrnd() % 8);

    // no hairy furballs for lofi version
    if (LOFI)
    {
      is_good = 0;
    }

    // randomises colour and position
    float r = (float)(xrnd() % 256) / 256.0f;
    float g = (float)(xrnd() % 256) / 256.0f;
    float b = (float)(xrnd() % 256) / 256.0f;
    float x = ((float)(xrnd() % 1024) / 1024.0f) * kWorldSize;
    float y = kHeight;
    float z = ((float)(xrnd() % 1024) / 1024.0f) * kWorldSize;

    // creates the instance
    ballz[i] = spawn_furball(x, y, z, is_good ? ultimate_furball : casual_furball, 8, is_good, r, g, b);

    // sets up ai
    ballz[i]->iq = 0;
    ballz[i]->smart = 0.05f + ((float)(xrnd() % 256) / 256.0f) * 0.15f;

    // adjusts ai for hairy furball
    if (is_good)
    {
      ballz[i]->smart *= 0.2f;
    }
    if (!is_good)
    {
      ballz[i]->scale = 0.2f + ((float)(xrnd() % 256) / 256.0f) * 0.15f;
    }
    else
    {
      ballz[i]->scale = 1.0f + ((float)(xrnd() % 256) / 256.0f) * 3.0f;
    }
    ballz[i]->dying = 0;
  }
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////  MECHANINCS FUNCTIONS
////////////////////////////////////////////////////

// creates [num] blood particles at [b] buffer
static void
explode(Blood* b, size_t num, float x, float y, float z)
{
  for (size_t i = 0; i < num; ++i)
  {
    // sets position
    b[i].x = x;
    b[i].y = y;
    b[i].z = z;
    // random velocity and colour
    b[i].vx = (((float)(xrnd() % 512) / 256.0f) - 1.0f) * 4.0f;
    b[i].vy = (((float)(xrnd() % 512) / 256.0f) - 1.0f) * 8.0f;
    b[i].vz = (((float)(xrnd() % 512) / 256.0f) - 1.0f) * 4.0f;
    b[i].r = 0.6f + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * 0.4f);
    b[i].g = 0.0f;
    b[i].b = 0.0f;
    b[i].alive = true;
  }
}

// shoots gun
static int
shoot(void)
{
  int32_t h = 0;
  Furball* hitball = NULL;
  float rhit[3] = {0.0f};

  play_sample(shot, 255, 128, 1000, 0);
  // gets shoot position and direction
  Vec3 m = {.x = cosf(looky) * cosf(lookx), .y = -sinf(lookx), .z = sinf(looky) * cosf(lookx)};
  float pdist = 3.40282347e+38F; // XXX: MAXFLOAT
  for (size_t i = 0; i < kBallz; ++i)
  {
    // picks furball closest to the ray
    if (ballz[i]->exists)
    {
      Vec3 p = {
        .x = ballz[i]->x,
        .y = ballz[i]->y,
        .z = ballz[i]->z};
      Vec3 p_minus_b = {
        .x = p.x - player.x,
        .y = p.y - player.y,
        .z = p.z - player.z};
      float t0 = vec3_dot(m, p_minus_b) / vec3_dot(m, m);
      if (t0 > 0)
      {
        rhit[0] = p.x - (player.x + t0 * m.x);
        rhit[1] = p.y - (player.y + t0 * m.y);
        rhit[2] = p.z - (player.z + t0 * m.z);
        if (length3v(rhit) < 10.0f)
        {
          if (vec3_dist(player, p) < pdist)
          {
            pdist = vec3_dist(player, p);
            h = 1;
            hitball = ballz[i];
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
    explode(hitball->blood, kParticles, hitball->x, hitball->y, hitball->z);

    // plays sound
    play_sample(slash, 255, 128, 900 + (xrnd() % 200), 0);
    play_sample(die[xrnd() % 8], 255, 128, 1000, 0);
  }
  return h;
}

// updates blood particles
static void
update_blood(Blood* b, size_t num)
{
  for (size_t i = 0; i < num; ++i)
  {
    // if the particle is alive...
    if (b[i].alive)
    {
      // updates position
      b[i].x += b[i].vx;
      b[i].y += b[i].vy;
      b[i].z += b[i].vz;

      // adds gravity
      b[i].vy -= 0.6f;

      // if it hits the ground
      if (b[i].y < 0.0f)
      {
        // no more updates for this one
        b[i].alive = false;

        // if was inside world sware, colours a pixel red
        int32_t pic_x = ((int32_t)b[i].x * ground_bmp->w) / (int32_t)kWorldSize;
        int32_t pic_y = ((int32_t)b[i].z * ground_bmp->h) / (int32_t)kWorldSize;
        if (pic_x >= 0 && pic_y >= 0 && pic_x < ground_bmp->w && pic_y < ground_bmp->h)
        {
          int32_t pix = getpixel(ground_bmp, pic_x, pic_y);
          int32_t pr = (int32_t)(((float)getr(pix) + b[i].r * 128.0f) / 1.5f);
          int32_t pg = (int32_t)(((float)getg(pix) + b[i].g * 128.0f) / 1.5f);
          int32_t pb = (int32_t)(((float)getb(pix) + b[i].b * 128.0f) / 1.5f);
          putpixel(ground_bmp, pic_x, pic_y, makecol(pr, pg, pb));
        }
      }
    }
  }
}

// updates furballs
static void
update_ballz(void)
{
  for (size_t i = 0; i < kBallz; ++i)
  {
    // if it's still alive
    if (ballz[i]->exists)
    {
      Furball* furball = ballz[i];

      // decreases direction change counter
      // and changes movement direction if needed
      furball->iq -= furball->smart;
      if (furball->iq <= 0)
      {
        furball->iq = (kIQ + (float)(xrnd() % kIQ)) * 3; // kPi;
        furball->speed = ((float)(xrnd() % 256) / 256.0f) * kSpeed;
        furball->a = ((float)(xrnd() % 256) / 256.0f) * kPiMultipliedBy2;
        furball->bounce_rate = 10.0f + ((float)(xrnd() % 256) / 256.0f) * 40.0f;
      }

      // adjusts bounce rate
      furball->bounce += 0.1f * kSpeed;

      // updates position
      furball->x += cosf(furball->a) * furball->speed;
      furball->z += sinf(furball->a) * furball->speed;
      furball->y = furball->scale * (furball->ultimate ? 1 : 8) + ABS(sinf(furball->iq) * furball->bounce_rate);

      // keeps inside the world
      if (furball->x < FLT_EPSILON)
      {
        furball->x = 0;
      }
      if (furball->z < FLT_EPSILON)
      {
        furball->z = 0;
      }
      if (furball->x > kWorldSize)
      {
        furball->x = kWorldSize;
      }
      if (furball->z > kWorldSize)
      {
        furball->z = kWorldSize;
      }

      // wites down trail for hair rendering
      for (size_t j = kTrail; j > 0; --j)
      {
        furball->trail[j][0] = furball->trail[j - 1][0];
        furball->trail[j][1] = furball->trail[j - 1][1];
        furball->trail[j][2] = furball->trail[j - 1][2];
      }
      furball->trail[0][0] = furball->x;
      furball->trail[0][1] = furball->y;
      furball->trail[0][2] = furball->z;
    }
    else if (ballz[i]->dying)
    {
      // if it's dying, only particles get updated
      update_blood(ballz[i]->blood, kParticles);
      ballz[i]->dying--;
    }
  }
}

// generates world map from 3 files
// first is colour, where it stores ground and grass colours
// second is grass height
// third is entity distribution, a specific colour represents one entity
// made so for quick level editing
static Buffer*
generate_world_map(const char* grass_file, const char* height_file, const char* entity_file, size_t density)
{
  Vec3 grass = {.x = 0.0f, .y = 0.0f, .z = 0.0f};

  // randomisation constants
  float colour_deviation = 0.1f;

  BITMAP* bmpg = NULL;
  BITMAP* bmph = NULL;
  // lofi gets smaller bitmaps
  if (LOFI)
  {
    BITMAP* bmpg_t = load_bmp(grass_file, 0);
    BITMAP* bmph_t = load_bmp(height_file, 0);
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
  int32_t size_x = bmpg->w;
  int32_t size_y = bmpg->h;

  // allocates grass buffer
  Buffer* b = malloc(sizeof(Buffer));
  b->vtx = malloc((size_t)(size_x * size_y) * density * 6 * sizeof(float));
  b->clr = malloc((size_t)(size_x * size_y) * density * 6 * sizeof(float));
  b->hardware = 0;
  size_t line_counter = 0;
  float* vertex = b->vtx;
  float* colour = b->clr;
  float patch_size = kWorldSize / (float)bmpg->w;

  for (int32_t x = 0; x < size_x; ++x)
  {
    for (int32_t y = 0; y < size_y; ++y)
    {
      for (size_t d = 0; d < density; ++d)
      {
        // for each bitmap pixel it creates [density] grass blades
        // colour and height get randomised a bit
        grass.x = (((float)x * kWorldSize) / (float)size_x) + ((((float)(xrnd() % 512) / 512.0f)) * patch_size);
        grass.y = (((float)y * kWorldSize) / (float)size_x) + ((((float)(xrnd() % 512) / 512.0f)) * patch_size);
        grass.z = (0.1f * (float)getr(getpixel(bmph, x, y))) + ((((float)(xrnd() % 512) / 512.0f)) * (grass.z * 0.3f));
        vertex[0] = grass.x;
        vertex[1] = 0;
        vertex[2] = grass.y;
        vertex[3] = grass.x;
        vertex[4] = grass.z;
        vertex[5] = grass.y;
        int32_t col = getpixel(bmpg, x, y);
        colour[3] = ((float)(getr(col)) / 256.0f) + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * colour_deviation);
        colour[4] = ((float)(getg(col)) / 256.0f) + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * colour_deviation);
        colour[5] = ((float)(getb(col)) / 256.0f) + ((((float)(xrnd() % 512) / 256.0f) - 1.0f) * colour_deviation);
        colour[0] = colour[3] * 0.5f;
        colour[1] = colour[4] * 0.5f;
        colour[2] = colour[5] * 0.5f;
        line_counter += 2;
        vertex += 6;
        colour += 6;
      }
    }
  }

  b->mode = GL_LINES;
  b->size = (int32_t)line_counter;
  vboize(b);
  ground_bmp = bmpg;

  destroy_bitmap(bmph);

  // creates entities
  bmph = load_bmp(entity_file, 0);
  size_x = bmph->w;
  size_y = bmph->h;
  for (int32_t x = 0; x < size_x; ++x)
  {
    for (int32_t y = 0; y < size_y; ++y)
    {
      // iterates through pixel looking for values representing entities
      grass.x = (((float)x * kWorldSize) / (float)size_x) + ((((float)(xrnd() % 512) / 512.0f)) * patch_size);
      grass.y = (((float)y * kWorldSize) / (float)size_x) + ((((float)(xrnd() % 512) / 512.0f)) * patch_size);
      grass.z = 1.0f + ((((float)(xrnd() % 512) / 512.0f)) * 1.0f);
      int32_t col = getpixel(bmph, x, y);
      if (col == makecol(0, 255, 0))
      {
        ents[num_ents].buf = tall_tree;
        ents[num_ents].x = grass.x;
        ents[num_ents].y = 0.0f;
        ents[num_ents].z = grass.y;
        ents[num_ents].s = grass.z * 4.0f;
        ents[num_ents].a = ((((float)(xrnd() % 512) / 512.0f)) * 360.0f);
        num_ents++;
      }
      if (col == makecol(255, 255, 128))
      {
        ents[num_ents].buf = cactus;
        ents[num_ents].x = grass.x;
        ents[num_ents].y = 0.0f;
        ents[num_ents].z = grass.y;
        ents[num_ents].s = grass.z * 2.0f;
        ents[num_ents].a = (((float)(xrnd() % 512) / 512.0f)) * 360.0f;
        num_ents++;
      }
      if (col == makecol(0, 128, 0))
      {
        ents[num_ents].buf = palm;
        ents[num_ents].x = grass.x;
        ents[num_ents].y = 0.0f;
        ents[num_ents].z = grass.y;
        ents[num_ents].s = grass.z * 3.0f;
        ents[num_ents].a = (((float)(xrnd() % 512) / 512.0f)) * 360.0f;
        num_ents++;
      }
      if (col == makecol(128, 64, 64))
      {
        ents[num_ents].buf = stick;
        ents[num_ents].x = grass.x;
        ents[num_ents].y = 0.0f;
        ents[num_ents].z = grass.y;
        ents[num_ents].s = grass.z * 1.5f;
        ents[num_ents].a = (((float)(xrnd() % 512) / 512.0f)) * 360.0f;
        num_ents++;
      }
      if (col == makecol(0, 255, 255))
      {
        ents[num_ents].buf = pine;
        ents[num_ents].x = grass.x;
        ents[num_ents].y = 0.0f;
        ents[num_ents].z = grass.y;
        ents[num_ents].s = grass.z * 3.0f;
        ents[num_ents].a = (((float)(xrnd() % 512) / 512.0f)) * 360.0f;
        num_ents++;
      }
      if (col == makecol(255, 0, 0))
      {
        ents[num_ents].buf = hut;
        ents[num_ents].x = grass.x;
        ents[num_ents].y = 0.0f;
        ents[num_ents].z = grass.y;
        ents[num_ents].s = grass.z * 2.0f;
        ents[num_ents].a = ((float)(xrnd() % 4)) * 90.0f;
        num_ents++;
      }
      if (col == makecol(255, 128, 128))
      {
        ents[num_ents].buf = church;
        ents[num_ents].x = grass.x;
        ents[num_ents].y = 0.0f;
        ents[num_ents].z = grass.y;
        ents[num_ents].s = grass.z * 6.0f;
        // ents[num_ents].a = 90.0f; //XXX: Original value
        ents[num_ents].a = (((float)(xrnd() % 512) / 512.0f)) * 360.0f;
        num_ents++;
      }
      if (col == makecol(128, 0, 0))
      {
        ents[num_ents].buf = brickhouse;
        ents[num_ents].x = grass.x;
        ents[num_ents].y = 0.0f;
        ents[num_ents].z = grass.y;
        ents[num_ents].s = grass.z * 1.5f;
        ents[num_ents].a = ((float)(xrnd() % 4)) * 90.0f;
        num_ents++;
      }
    }
  }
  destroy_bitmap(bmph);
  return b;
}

// generate everything batch
static void
generate_stuff(void)
{
  // generate_font("cartographer.ttf",1000);
  ultimate_furball = generate_furball_ultimate(1.0f, kTrail, 8);
  casual_furball = generate_furball_normal(8.0f, kTrail, 8);
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
  g_grass = generate_world_map("grassmap.bmp", "grassheight.bmp", "entitymap.bmp", LOFI ? 2 : 4);
  create_ballz();
}

// timed function that updates everything
// and handles input
static void
timer_proc(void)
{
  float move_spd = kMoveSpeed;

  // if it's intro, waits for space
  if (game_state == INTRO)
  {
    position_mouse(kScreenWidth / 2, kScreenHeight / 2);
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
        talk_counter = kTalkDelay + (xrnd() % kTalkDelay);
        // play_sample(furtalk[(xrnd()%16)>>1],128,128,1000,0);
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
        for (size_t i = 0; i < kBallz; ++i)
        {
          ballz[i]->a = atan2f(-ballz[i]->z + player.z, -ballz[i]->x + player.x);
          ballz[i]->iq += kPiMultipliedBy4;
        }
      }
    }
    // shift makes you run faster
    if (key[KEY_LSHIFT])
    {
      move_spd *= 2;
    }

    if (key[KEY_ESC])
    {
      game_over = true;
    }

    // rotates camera using mouse input
    get_mouse_mickeys(&mx, &my);
    position_mouse(kScreenWidth / 2, mouse_y);
    update_blood(g_blood, kParticles);
    lookx = kPiDividedBy2 * (float)(mouse_y - (kScreenHeight / 2)) / (float)(kScreenHeight / 2);
    looky += kPi * (float)mx / (kScreenWidth / 2.0f);

    // moves player according to key input
    if (key[KEY_LEFT] || key[KEY_A])
    {
      player.x += sinf(looky) * move_spd;
      player.z += -cosf(looky) * move_spd;
      state = WALK;
    }
    if (key[KEY_RIGHT] || key[KEY_D])
    {
      player.x += -sinf(looky) * move_spd;
      player.z += cosf(looky) * move_spd;
      state = WALK;
    }

    if (key[KEY_UP] || key[KEY_W])
    {
      player.x += cosf(looky) * move_spd;
      player.z += sinf(looky) * move_spd;
      state = WALK;
    }
    if (key[KEY_DOWN] || key[KEY_S])
    {
      player.x -= cosf(looky) * move_spd;
      player.z -= sinf(looky) * move_spd;
      state = WALK;
    }

    // if not walking, sets idle state
    if (!key[KEY_UP] && !key[KEY_DOWN] && state != IDLE)
    {
      state = IDLE;
    }

    // if not moving, add breathing effect
    if (state == STOP)
    {
      if (ABS(player.y - kHeight) > kHaystack * sinf((player.x / kWorldSize) * kPi) * sinf((player.x / kWorldSize) * kPi) * sinf((player.z / kWorldSize) * kPi) * sinf((player.z / kWorldSize) * kPi))
      {
        player.y -= SGN(player.y - kHeight) * 0.4f;
      }
      else
      {
        state = IDLE;
        lookf = 0.0f;
      }
    }
    // shoots
    if (mouse_b)
    {
      shoot();
    }
  }

  bounce = kBounceRate;

  // if idle, add breathe effect aswell
  if (state == IDLE)
  {
    lookf += 0.2f;
    player.y = kHeight + (bounce * 0.3f * (sinf(lookf)));
  }
  if (state == WALK)
  {
    lookf += 0.1f;

    player.y = kHeight + (bounce * (sinf(lookf)));
  }

  // keeps player inside the world
  if (player.x > (kWorldSize - 50.0f))
  {
    player.x = (kWorldSize - 50.0f);
  }
  if (player.z > (kWorldSize - 50.0f))
  {
    player.z = (kWorldSize - 50.0f);
  }
  if (player.x < 50.0f)
  {
    player.x = 50.0f;
  }
  if (player.z < 50.0f)
  {
    player.z = 50.0f;
  }

  // updates furballs
  update_ballz();
}

// rendering function
static void
draw(void)
{
  char nums[16] = {0};
  int32_t left = 0;

  // checks for number of furballs left
  for (size_t i = 0; i < kBallz; ++i)
  {
    if (ballz[i]->exists)
    {
      left++;
    }
  }

  // a cheatcode ;)
  if (key[KEY_Y])
  {
    left = 0;
  }

  // sets state to outro if all furballs are dead
  if (!left)
  {
    game_state = OUTRO;
  }

  // stop furball talking if you killed too much
  if (left < 100)
  {
    talking = 0;
  }
  sprintf(nums, "%03i:%03i", left, kBallz);

  // up sets furball count display
  for (size_t i = 0; i < 7; ++i)
  {
    nums[i] -= 0x30;
  }

  // basic gl frame setup
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glPushMatrix();
  glFrustum(-(float)(kScreenWidth / kScreenHeight), (float)(kScreenWidth / kScreenHeight), -1.0f, 1.0f, 1.0f, kDrawDist);
  glMatrixMode(GL_MODELVIEW);

  // Clear the RGB buffer and the depth buffer
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Set the modelview matrix to be the identity matrix
  glLoadIdentity();
  glPushMatrix();

  // Set the camera
  glRotatef(DEG(lookx), 1.0f, 0.0f, 0.0f);
  glRotatef(DEG(looky) + 90.0f, 0.0f, 1.0f, 0.0f);
  glRotatef(DEG(lookz), 0.0f, 0.0f, 1.0f);
  glTranslatef(-player.x, -player.y, -player.z);

  // Save the camera matrix
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glTranslatef(0.0f, sinf(lookx * 1.333f), 0.0f);
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glDisable(GL_DEPTH_TEST);
  glDepthMask(GL_FALSE);

  // draws backgorund
  glBegin(GL_QUADS);
  glColor3f(1.0f, 0.8f, 0.4f);
  glVertex2f(-1.0f, 0.5f);
  glVertex2f(1.0f, 0.5f);
  glColor3f(1.0f, 0.85f, 0.6f);
  glVertex2f(1.0f, 0.2f);
  glVertex2f(-1.0f, 0.2f);

  glVertex2f(-1.0f, 0.2f);
  glVertex2f(1.0f, 0.2f);
  glColor3f(1.0f, 0.9f, 0.8f);
  glVertex2f(1.0f, 0.0f);
  glVertex2f(-1.0f, 0.0f);

  glVertex2f(-1.0f, 0.0f);
  glVertex2f(1.0f, 0.0f);
  glColor3f(0.0f, 0.0f, 0.0f);
  glVertex2f(1.0f, -1.0f);
  glVertex2f(-1.0f, -1.0f);
  glEnd();

  glDepthMask(GL_TRUE);
  glEnable(GL_DEPTH_TEST);
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();

  // Translate and rotate the object
  glColor3f(1.0f, 1.0f, 1.0f);

  // draws world
  draw_tree();

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();

  glDisable(GL_TEXTURE_2D);
  glColor3f(1.0f, 1.0f, 1.0f);
  mx = (mouse_x - 160) / 160;
  my = -(mouse_y - 120) / 120;

  glPointSize(3.0f);
  glDisable(GL_DEPTH_TEST);
  if (game_state == PLAY)
  {

    // draws crosshair
    glBegin(GL_POINTS);
    glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
    glVertex2f(0.0f, 0.0f);
    glEnd();

    // draws numbers
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, numbers);
    glBegin(GL_QUADS);

    float ns = 0.0f;
    for (size_t i = 0; i < 7; ++i, ns += kNumSize)
    {

      glTexCoord2f(nums[i] * 0.090909f, 1.0f);
      glVertex2f(-1.0f + ns, 1.0f);
      glTexCoord2f((float)(nums[i] + 1) * 0.090909f, 1.0f);
      glVertex2f(-1.0f + ns + (kNumSize * 1.2f), 1.0f);
      glTexCoord2f((float)(nums[i] + 1) * 0.090909f, 0.0f);
      glVertex2f(-1.0f + ns + (kNumSize * 1.2f), 1.0f - kNumSize);
      glTexCoord2f(nums[i] * 0.090909f, 0.0f);
      glVertex2f(-1.0f + ns, 1.0f - kNumSize);
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
    glTexCoord2f(0.0f, 1.0f);
    glVertex2f(-1.0f, 1.0f);
    glTexCoord2f(1.0f, 1.0f);
    glVertex2f(1.0f, 1.0f);
    glTexCoord2f(1.0f, 0.0f);
    glVertex2f(1.0f, -1.0f);
    glTexCoord2f(0.0f, 0.0f);
    glVertex2f(-1.0f, -1.0f);
    glEnd();
    glDisable(GL_TEXTURE_2D);
  }

  // renders outro screen
  if (game_state == OUTRO)
  {
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, outro);
    glBegin(GL_QUADS);
    glTexCoord2f(0.0f, 1.0f);
    glVertex2f(-1.0f, 1.0f);
    glTexCoord2f(1.0f, 1.0f);
    glVertex2f(1.0f, 1.0f);
    glTexCoord2f(1.0f, 0.0f);
    glVertex2f(1.0f, -1.0f);
    glTexCoord2f(0.0f, 0.0f);
    glVertex2f(-1.0f, -1.0f);
    glEnd();
    glDisable(GL_TEXTURE_2D);
  }
  glEnable(GL_DEPTH_TEST);

  // flip buffers
  allegro_gl_flip();
  // printf("Draw took %u msec\n",clock()-clk);
}

// timer thread
static void
timer(void)
{
  tim++;
}
END_OF_FUNCTION(timer)

// application entry point
int
main(int argc, char** argv)
{
  (void)argc;

  float fogc[] = {1.0f, 0.8f, 0.4f, 0.0f};
  char samp[64] = {0};

  //  allegro setup
  allegro_init();
  install_allegro_gl();

  allegro_gl_clear_settings();
  allegro_gl_set(AGL_COLOR_DEPTH, 32);
  allegro_gl_set(AGL_Z_DEPTH, 24);
  allegro_gl_set(AGL_WINDOWED, FALSE);
  allegro_gl_set(AGL_DOUBLEBUFFER, 1);
  allegro_gl_set(AGL_SUGGEST, AGL_COLOR_DEPTH | AGL_Z_DEPTH | AGL_DOUBLEBUFFER | AGL_FULLSCREEN);

  if (set_gfx_mode(GFX_OPENGL, kScreenWidth, kScreenHeight, 0, 0) < 0)
  {
    allegro_message("Error setting OpenGL graphics mode:\n%s\n"
                    "Allegro GL error : %s\n",
                    allegro_error,
                    allegro_gl_error);
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

  glFogf(GL_FOG_START, kDrawDist * (LOFI ? 0.8f : 0.6f));
  glFogf(GL_FOG_END, kDrawDist);
  glFogfv(GL_FOG_COLOR, (float*)&fogc);

  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);

  // timer setup
  install_int_ex(timer, BPS_TO_TIMER(30));
  // srand((unsigned int)time(NULL));
  // rngstate[0] = (uint64_t)time(NULL);
  xrnd_seed(time(NULL));

  // generates everything
  glEnable(GL_TEXTURE_2D);
  generate_stuff();

  // more gl setup
  glClearColor(1.0f, 0.8f, 0.4f, 0.0f);
  glEnable(GL_BLEND);
  glEnable(GL_ALPHA_TEST);
  glAlphaFunc(GL_GREATER, 0.3f);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  player.x = 20.0f;
  player.z = 20.0f;
  looky = kPiDividedBy4;

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

  for (size_t i = 0; i < kSamps; ++i)
  {
    sprintf(samp, "fur%zu.wav", i + 1);
    furtalk[i] = load_wav(samp);
    sprintf(samp, "die%zu.wav", i + 1);
    die[i] = load_wav(samp);
  }

  // sets player initial position and direction
  player.x = 150;
  player.y = kHeight;
  player.z = 150;
  lookx = kPiDividedBy2 * (float)((kScreenHeight / 2) - (kScreenWidth / 2)) / (float)(kScreenWidth / 2);
  looky = 1.639519f;
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
  } while (!game_over);

  set_gfx_mode(GFX_TEXT, 0, 0, 0, 0);

  return 0;
}
END_OF_MAIN()
