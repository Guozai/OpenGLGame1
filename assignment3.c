/*
 * This program shows how to do collision detection between
 * particles in two dimensions.
 *
 * Uses hard-sphere assumptions for collision dynamics.
 *
 * Both brute-force and uniform grid approaches implemented.
 *
 * Numerical integration of equations of motion.
 *
 * $Id: particles_2D.C,v 1.13 2012/09/24 00:57:22 gl Exp gl $
 *
 * Make sure you add the SDL2 library to the library list when compiling: -lSDL2
 */
#define GL_GLEXT_PROTOTYPES

#include <math.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>
#include <SDL2/SDL.h>

#if _WIN32
#  include <Windows.h>
#endif
#if __APPLE__
#  include <OpenGL/gl.h>
#  include <OpenGL/glu.h>
#else
#  include <GL/gl.h>
#  include <GL/glu.h>
#endif

#include "shaders.h"

#define ZOOM_SCALE 4.0
#define MAXV (-3.0) // the speed of absorb for smooth absorption
#define MAX_GAME_LEVEL 3

SDL_Window *mainWindow = 0;

/* Debugging controls .*/
enum debugFlags {
    debug_time,
    debug_wall,
    debug_initialise_particle,
    debug_particle,
    debug_particle_collision,
    debug_collideParticlesBruteForce,
    debug_collideParticlesUniformGrid,
    debug_collisionReactionParticles2DbasisChange,
    debug_collisionReactionParticles2DprojNormal,
    debug_framerate,
    debug_range_check,
    debug_sum_kinetic_energy,
    debug_sum_momentum,
    debug_inelastic,
    debug_mouse,
    numDebugFlags
};
int debug[numDebugFlags] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

/* Use our type so we can change precision easily. */
typedef double Real;

/* Small number to handle numerical imprecision. */
const Real epsilon = 1.0e-6;

typedef enum { inactive, zoom } CameraControl;

struct camera_t {
    int lastX, lastY;
    float scale;
    CameraControl control;
} camera = { 0, 0, ZOOM_SCALE, inactive };

/* Particles (particles). */
struct Particle {
    Real position[2];
    Real velocity[2];
    Real velocityF[2]; // for smooth absorption
    Real radius;
    Real radiusF; // final radius after smooth absorption
    Real t; // time the aborted sphere gets totally inside of the larger absorber
    Real count_h; // count the total of h bigger than particle.t or not
    Real vRad; // velocity of radius changing from radius to radiusF(inal)
    Real mass;
    Real elasticity;
    GLUquadric *quadric;  /* For rendering. */
    int slices, loops;    /* For rendering. */
    int isAbort; /* abort flag for inelastic collision */
};

/* Control 1D or 2D */
const int dimension = 2;

/* Control random or explicit initial positions */
enum {
    randomly,
    explicitly
} const initialiseParticles = randomly;

/* Control collision reaction calculation */
enum ReactionCalculation {
    basisChange,
    projNormal
} reacCalc = basisChange;

const int numParticles = 100;
const int maxParticles = 200;
Particle particle[maxParticles];

/* Arena. */
struct Arena {
    Real min[2], max[2];
    Real momentum[2];
} arena;

/* Rendering info. */
enum renderMode { wire, solid };
static const int milli = 1000;

// global parameters
struct global {
    int gameLevel;
    renderMode renMode;
    Real elapsedTime, startTime, pauseTime, resumeTime, pauseInterval;
    bool go, pause;
    bool eject;
    float x, y;
    Real maxRadius;
    int numMaxRadius;
    bool win, lose;
} g = { 1, wire, 0.0, 0.0, 0.0, 0.0, 0.0, false, false, false, 0.0, 0.0, 0.0, 0, false, false };

/* Collision detection method. */
enum CollisionDetectionMethod {
    bruteForce,
    uniformGrid
};

CollisionDetectionMethod CDmethod = bruteForce;

Real random_uniform() {
    return rand()/(float)RAND_MAX;
}

// uint variable for shader
GLuint shaderProgram;

int wantRedisplay = 1;
void postRedisplay()
{
    wantRedisplay = 1;
}

void panic(const char *m) {
    printf("%s", m);
    exit(1);
}

void initialiseArena()
{
    const Real halfLength = 7.5;
    
    arena.min[0] = -halfLength;
    arena.min[1] = -halfLength;
    arena.max[0] = halfLength;
    arena.max[1] = halfLength;
    
    arena.momentum[0] = 0.0;
    arena.momentum[1] = 0.0;
}

/* player */
void initialiseParticlesExplicitly()
{
    GLUquadric *quadric = gluNewQuadric();
    
    particle[0].position[0] = 0.0 * arena.max[0];
    particle[0].position[1] = 0.0 * arena.max[1];
    particle[0].velocity[0] = particle[0].velocityF[0] = 0.0;
    particle[0].velocity[1] = particle[0].velocityF[1] = 0.0;
    particle[0].mass = 0.05;
    particle[0].radius = particle[0].radiusF = sqrt(particle[0].mass);
    particle[0].t = 0.0;
    particle[0].count_h = 0.0;
    particle[0].vRad = 0.0;
    particle[0].elasticity = 1.0;
    particle[0].quadric = quadric;
    particle[0].slices = 10;
    particle[0].loops = 3;
    particle[0].isAbort = false;
    
    for (int i = numParticles; i < maxParticles; i++) {
        particle[i].position[0] = 1.0 * arena.max[0];
        particle[i].position[1] = 1.0 * arena.max[1];
        particle[i].velocity[0] = particle[i].velocityF[0] = 0.0;
        particle[i].velocity[1] = particle[i].velocityF[1] = 0.0;
        particle[i].mass = 0.005;
        particle[i].radius = particle[i].radiusF = sqrt(particle[0].mass);
        particle[i].t = 0.0;
        particle[i].count_h = 0.0;
        particle[i].vRad = 0.0;
        particle[i].elasticity = 1.0;
        particle[i].quadric = quadric;
        particle[i].slices = 10;
        particle[i].loops = 3;
        particle[i].isAbort = true;
    }
}

void initialiseParticlesRandomly(Real maxV)
{
    GLUquadric *quadric = gluNewQuadric();
    const Real maxVelocity = maxV;
    Real n[2], n_mag_sq, sum_radii, sum_radii_sq;
    bool collision, done;
    int i, j;
    
    for (i = 1; i < numParticles; i++) {
        particle[i].velocity[0] = particle[i].velocityF[0] = (random_uniform() - 0.5) * maxVelocity;
        particle[i].velocity[1] = particle[i].velocityF[1] = (random_uniform() - 0.5) * maxVelocity;
        particle[i].mass = random_uniform() * 0.07;
        
        particle[i].radius = particle[i].radiusF = sqrt(particle[i].mass);
        particle[i].t = 0.0;
        particle[i].count_h = 0.0;
        particle[i].vRad = 0.0;
        particle[i].elasticity = 1.0;
        particle[i].quadric = quadric;
        particle[i].slices = 10;
        particle[i].loops = 3;
        particle[i].isAbort = false;
        
        done = false;
        while (!done) {
            particle[i].position[0] = random_uniform() *
            (arena.max[0] - arena.min[0] - 2.0 * particle[i].radius) +
            arena.min[0] + particle[i].radius + epsilon;
            particle[i].position[1] = random_uniform() *
            (arena.max[1] - arena.min[1] - 2.0 * particle[i].radius) +
            arena.min[1] + particle[i].radius + epsilon;
            
            /* Check for collision with existing particles. */
            collision = false;
            j = 0;
            while (!collision && j < i) {
                sum_radii = particle[i].radius + particle[j].radius;
                sum_radii_sq = sum_radii * sum_radii;
                n[0] = particle[j].position[0] - particle[i].position[0];
                n[1] = particle[j].position[1] - particle[i].position[1];
                n_mag_sq = n[0] * n[0] + n[1] * n[1];
                if (n_mag_sq < sum_radii_sq)
                collision = true;
                else
                j++;
            }
            if (!collision)
            done = true;
        }
        if (debug[debug_initialise_particle])
        printf ("initialiseParticles: x %f y %f\n",
                particle[i].position[0], particle[i].position[1]);
    }
}

void setRenderMode(renderMode rm)
{
    /* Example of GNU C/C++ brace indentation style.  */
    if (rm == wire)
    {
        glDisable(GL_LIGHTING);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_NORMALIZE);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    }
    else if (rm == solid)
    {
        glEnable(GL_LIGHTING);
        glEnable(GL_LIGHT0);
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_NORMALIZE);
        glShadeModel(GL_SMOOTH);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        
    }
}

void changeRenderMode(void)
{
    if (g.renMode == wire) {
        g.renMode = solid;
    } else {
        g.renMode = wire;
    }
    setRenderMode(g.renMode);
}

void init()
{
    setRenderMode(g.renMode);
    initialiseArena();
    initialiseParticlesExplicitly(); //initial player mote
    initialiseParticlesRandomly(0.0);
    
    shaderProgram = getShader("shader.vert", "shader.frag");
}

void initStage(int stage)
{
    switch(stage) {
        case 2:
        initialiseParticlesExplicitly();
        initialiseParticlesRandomly(1.0);
        break;
        case 3:
        initialiseParticlesExplicitly();
        initialiseParticlesRandomly(3.0);
        break;
        case 1:
        default:
        initialiseParticlesExplicitly();
        initialiseParticlesRandomly(0.0);
        break;
    }
}

void displayArena()
{
    glBegin(GL_LINE_LOOP);
    glVertex3f(arena.min[0], arena.min[1], 0.0);
    glVertex3f(arena.max[0], arena.min[1], 0.0);
    glVertex3f(arena.max[0], arena.max[1], 0.0);
    glVertex3f(arena.min[0], arena.max[1], 0.0);
    glEnd();
}

void displayParticle(Particle *p, float sx, float sy, float sz)
{
    glPushMatrix();
    glScalef(sx, sy, sz);
    gluDisk(p->quadric, 0.0, p->radius, p->slices, p->loops);
    glPopMatrix();
}

float sumKineticEnergy()
{
    Real v_sq, K;
    
    K = 0;
    for (int i = 0; i < maxParticles; i++) {
        v_sq = particle[i].velocity[0] * particle[i].velocity[0] +
        particle[i].velocity[1] * particle[i].velocity[1];
        K += 0.5 * particle[i].mass * v_sq;
    }
    
    return K;
}

void sumMomentum(Real *p)
{
    p[0] = p[1] = 0;
    for (int i = 0; i < maxParticles; i++) {
        p[0] += particle[i].mass * particle[i].velocity[0];
        p[1] += particle[i].mass * particle[i].velocity[1];
    }
    p[0] += arena.momentum[0];
    p[1] += arena.momentum[1];
}

void eulerStepSingleParticle(Particle &p, Real h)
{
    p.position[0] += h * p.velocity[0];
    p.position[1] += h * p.velocity[1];
}

void eulerStepSmoothAbsorption(Particle &p, Real h)
{
    p.position[0] += h * p.velocity[0];
    p.position[1] += h * p.velocity[1];
    if (p.radiusF > p.radius) {
        p.radius += h * p.vRad;
    }
    if (p.isAbort) {
        p.count_h += h;
        if (p.isAbort && p.count_h > p.t) {
            p.velocity[0] = p.velocity[1] = 0.0;
            p.radius = p.radiusF = 0.0;
            p.count_h = 0.0;
        }
    }
}

void integrateMotionParticles(Real h)
{
    for (int i = 0; i < maxParticles; i++)
    eulerStepSmoothAbsorption(particle[i], h);
}

void collideParticleWall(Particle &p, Arena &a)
{
    float dp[2];
    
    dp[0] = dp[1] = 0.0;
    if ((p.position[0] - p.radius) < a.min[0]) {
        p.position[0] +=
        2.0 * (a.min[0] - (p.position[0] - p.radius));
        p.velocity[0] *= -1.0;
        dp[0] += p.mass * -2.0 * p.velocity[0];
    }
    if ((p.position[1] - p.radius) < a.min[1]) {
        p.position[1] +=
        2.0 * (a.min[1] - (p.position[1] - p.radius));
        p.velocity[1] *= -1.0;
        dp[1] += p.mass * -2.0 * p.velocity[1];
    }
    if ((p.position[0] + p.radius) > a.max[0]) {
        p.position[0] -=
        2.0 * (p.position[0] + p.radius - a.max[0]);
        p.velocity[0] *= -1.0;
        dp[0] += p.mass * -2.0 * p.velocity[0];
    }
    if ((p.position[1] + p.radius) > a.max[1]) {
        p.position[1] -=
        2.0 * (p.position[1] + p.radius - a.max[1]);
        p.velocity[1] *= -1.0;
        dp[1] += p.mass * -2.0 * p.velocity[1];
    }
    arena.momentum[0] += dp[0];
    arena.momentum[1] += dp[1];
}

void collideParticlesWall()
{
    for (int i = 0; i < maxParticles; i++) {
        if (debug[debug_wall])
        printf("%d %f %f\n",
               i, particle[i].position[0], particle[i].position[1]);
        collideParticleWall(particle[i], arena);
    }
}

inline bool collisionDetectionParticles(Particle &p1, Particle &p2)
{
    Real sum_radii, sum_radii_sq, n[2], n_mag_sq;
    
    sum_radii = p1.radius + p2.radius;
    sum_radii_sq = sum_radii * sum_radii;
    n[0] = p2.position[0] - p1.position[0];
    n[1] = p2.position[1] - p1.position[1];
    n_mag_sq = n[0] * n[0] + n[1] * n[1];
    if (!p1.isAbort && !p2.isAbort && n_mag_sq <= sum_radii_sq)
    return true;
    else
    return false;
}

void collisionReactionParticles1D(Particle &p1, Particle &p2)
{
    float m1, m2, v1i, v2i, vf;
    
    m1 = p1.mass;
    m2 = p2.mass;
    v1i = p1.velocity[0];
    v2i = p2.velocity[0];
    vf = (m1 * v1i + m2 * v2i) / (m1 + m2);
    
    p1.radiusF = cbrt((m1 + m2) / m2) * p1.radius;
    // for smooth absorption
    if (p1.radius > p2.radius)
    p2.t = 2.0 * p2.radius / -MAXV;
    else
    p2.t = 2.0 * p1.radius / -MAXV;
    p1.vRad = (p1.radiusF - p1.radius) / p2.t;
    p1.mass = m1 + m2;
    p1.velocityF[0] = vf;
    p2.velocity[0] = MAXV;
    p2.velocity[1] = 0.0;
    if (debug[debug_inelastic])
    p2.radius = 0.1;
    p2.isAbort = true;
}

void collisionReactionParticles2DprojNormal(Particle &p1, Particle &p2)
{
    Real n[2], n_mag;
    Real projnv1, projnv2;
    Real m1, m2, v1i, v2i, vf;
    
    /* Normal vector n between centres. */
    n[0] = p2.position[0] - p1.position[0];
    n[1] = p2.position[1] - p1.position[1];
    
    /* Normalise n. */
    n_mag = sqrt(n[0] * n[0] + n[1] * n[1]);
    n[0] /= n_mag;
    n[1] /= n_mag;
    
    /* Vector projection/component/resolute of velocity in n direction. */
    projnv1 = n[0] * p1.velocity[0] + n[1] * p1.velocity[1];
    projnv2 = n[0] * p2.velocity[0] + n[1] * p2.velocity[1];
    
    /* Use 1D equations to calculate final velocities in n direction. */
    v1i = projnv1;
    v2i = projnv2;
    m1 = p1.mass;
    m2 = p2.mass;
    vf = (m1 * v1i + m2 * v2i) / (m1 + m2);
    
    /* Vector addition to solve for final velocity. */
    p1.radiusF = cbrt((m1 + m2) / m1) * p1.radius; // combined mote to be of size proportional to total mass
    // for smooth absorption
    if (p1.radius > p2.radius)
    p2.t = 2.0 * p2.radius / -MAXV;
    else
    p2.t = 2.0 * p1.radius / -MAXV;
    p1.vRad = (p1.radiusF - p1.radius) / p2.t;
    p1.velocityF[0] = vf * n[0];
    p1.velocityF[1] = vf * n[1];
    p1.mass = m1 + m2;
    Real theta = atan2f((p2.position[1] - p1.position[1]), (p2.position[0] - p1.position[0]));
    p2.velocity[0] = MAXV * cosf(theta);
    p2.velocity[1] = MAXV * sinf(theta);
    if (debug[debug_inelastic])
    p2.radius = 0.1;
    p2.isAbort = true;
}

void collisionReactionParticles2DbasisChange(Particle &p1, Particle &p2)
{
    Real n[2], t[2], n_mag;
    Real v1_nt[2], v2_nt[2];
    Real m1, m2, v1i, v2i, vf;
    
    if (debug[debug_collisionReactionParticles2DbasisChange]) {
        printf("collisionReactionParticles2DbasisChange:\n");
        printf("velocities before: %f %f %f %f\n",
               p1.velocity[0], p1.velocity[1],
               p2.velocity[0], p2.velocity[1]);
    }
    
    /* Normal vector n between centres. */
    n[0] = p2.position[0] - p1.position[0];
    n[1] = p2.position[1] - p1.position[1];
    
    /* Normalise n. */
    n_mag = sqrt(n[0] * n[0] + n[1] * n[1]);
    n[0] /= n_mag;
    n[1] /= n_mag;
    
    /* Tangent vector t. */
    t[0] = -n[1];
    t[1] = n[0];
    
    /* Change basis for velocities from standard basis to nt basis. */
    v1_nt[0] = n[0] * p1.velocity[0] + n[1] * p1.velocity[1];
    v1_nt[1] = t[0] * p1.velocity[0] + t[1] * p1.velocity[1];
    v2_nt[0] = n[0] * p2.velocity[0] + n[1] * p2.velocity[1];
    v2_nt[1] = t[0] * p2.velocity[0] + t[1] * p2.velocity[1];
    
    /* Use 1D equations to calculate final velocities in n direction. */
    m1 = p1.mass;
    m2 = p2.mass;
    v1i = v1_nt[0];
    v2i = v2_nt[0];
    vf = (m1 * v1i + m2 * v2i) / (m1 + m2);
    
    /* Update the 2D velocity. Force is in n direction, so in nt basis,
     * velocity change is in n direction only, no change in t direction.
     */
    p1.radiusF = cbrt((m1 + m2) / m1) * p1.radius; // combined mote to be of size proportional to total mass
    // for smooth absorption
    if (p1.radius > p2.radius)
    p2.t = 2.0 * p2.radius / -MAXV;
    else
    p2.t = 2.0 * p1.radius / -MAXV;
    p1.vRad = (p1.radiusF - p1.radius) / p2.t;
    p1.mass = m1 + m2;
    v1_nt[0] = vf;
    Real theta = atan2f((p2.position[1] - p1.position[1]), (p2.position[0] - p1.position[0]));
    p2.velocity[0] = MAXV * cosf(theta);
    p2.velocity[1] = MAXV * sinf(theta);
    if (debug[debug_inelastic])
    p2.radius = 0.1;
    p2.isAbort = true;
    
    /* Change back to standard basis. */
    p1.velocityF[0] = n[0] * v1_nt[0] + t[0] * v1_nt[1];
    p1.velocityF[1] = n[1] * v1_nt[0] + t[1] * v1_nt[1];
    // p2.velocity[0] = n[0] * v2_nt[0] + t[0] * v2_nt[1];
    // p2.velocity[1] = n[1] * v2_nt[0] + t[1] * v2_nt[1];
    
    if (debug[debug_collisionReactionParticles2DbasisChange]) {
        printf("velocities after: %f %f %f %f\n",
               p1.velocity[0], p1.velocity[1],
               p2.velocity[0], p2.velocity[1]);
    }
}

void collideParticlesBruteForce(Real h)
{
    int i, j;
    
    for (i = 0; i < maxParticles - 1; i++) {
        if (!particle[i].isAbort) {
            for (j = i + 1; j < maxParticles; j++) {
                if (!particle[j].isAbort) {
                    if (collisionDetectionParticles(particle[i], particle[j])) {
                        if (g.gameLevel == 1 && i == 0 && particle[i].mass < particle[j].mass) {
                            g.lose = true;
                            g.go = false;
                            printf("Player loses!\nPress [SPACE] bar to Restart the level...\n");
                        }
                        
                        if (debug[debug_collideParticlesBruteForce])
                        printf("collideParticlesBruteForce: collision %d %d\n", i, j);
                        
                        /* Take step back. Better approaches possible. */
                        // eulerStepSingleParticle(particle[i], -h);
                        // eulerStepSingleParticle(particle[j], -h);
                        
                        if (debug[debug_collideParticlesBruteForce]) {
                            // printf("velocities before: %f %f %f %f\n",
                            //        particle[i].velocity[0], particle[i].velocity[1],
                            //        particle[j].velocity[0], particle[j].velocity[1]);
                            printf("mass before: %f %f\n",
                                   particle[i].mass, particle[j].mass);
                        }
                        
                        /* Collision */
                        if (debug[debug_inelastic]){
                            printf("Bf particle %d: radius %f, mass %f  particle %d: radius %f, mass %f\n", i, particle[i].radius, particle[i].mass, j, particle[j].radius, particle[j].mass);
                        }
                        if (dimension == 1) {
                            collisionReactionParticles1D(particle[i], particle[j]);
                        }
                        else if (dimension == 2) {
                            if (reacCalc == basisChange)
                            collisionReactionParticles2DbasisChange(particle[i], particle[j]);
                            else if (reacCalc == projNormal)
                            collisionReactionParticles2DprojNormal(particle[i], particle[j]);
                            else
                            panic("collision reaction calculation not specified\n");
                        }
                        if (debug[debug_inelastic])
                        printf("Af particle %d: radius %f, mass %f  particle %d: radius %f, mass %f\n", i, particle[i].radius, particle[i].mass, j, particle[j].radius, particle[j].mass);
                        
                        if (debug[debug_collideParticlesBruteForce]) {
                            // printf("velocities after: %f %f %f %f\n",
                            //        particle[i].velocity[0], particle[i].velocity[1],
                            //        particle[j].velocity[0], particle[j].velocity[1]);
                            printf("mass after: %f %f\n",
                                   particle[i].mass, particle[j].mass);
                        }
                        
                        /* Step forward. */
                        // eulerStepSingleParticle(particle[i], h);
                        // eulerStepSingleParticle(particle[j], h);
                        
                        // if (debug[debug_collideParticlesBruteForce]) {
                        // printf("velocities after: %f %f %f %f\n",
                        //        particle[i].velocity[0], particle[i].velocity[1],
                        //        particle[j].velocity[0], particle[j].velocity[1]);
                        // }
                        
                        /* Check walls. */
                        collideParticleWall(particle[i], arena);
                        collideParticleWall(particle[j], arena);
                        
                        // if (debug[debug_collideParticlesBruteForce]) {
                        // printf("velocities after: %f %f %f %f\n",
                        //        particle[i].velocity[0], particle[i].velocity[1],
                        //        particle[j].velocity[0], particle[j].velocity[1]);
                        // }
                    }
                    if (g.maxRadius < particle[j].radius)
                    g.maxRadius = particle[j].radius;
                    if (g.maxRadius == particle[j].radius)
                    g.numMaxRadius++;
                }
            }
            if (g.maxRadius < particle[i].radius)
            g.maxRadius = particle[i].radius;
            if (g.maxRadius == particle[i].radius)
            g.numMaxRadius++;
        }
    }
}

inline void calcGridIndex(Particle &p, Arena a,
                          Real gridCellSize[2], int gridNumCells[2],
                          int index[2])
{
    index[0] = (int)((p.position[0] - a.min[0]) / gridCellSize[0]);
    index[1] = (int)((p.position[1] - a.min[1]) / gridCellSize[1]);
    if (debug_range_check == 0) {
        if (index[0] < 0 || index[0] > gridNumCells[0] - 1)
        panic("gridIndex: index out of range\n");
        if (index[1] < 0 || index[1] > gridNumCells[1] - 1)
        panic("gridIndex: index out of range\n");
    }
}

void collideParticlesUniformGrid(Real h)
{
    Real gridCellSize[2];
    int **gridCellParticleCount, **gridCellParticleListEnd, *gridCellParticleList;
    int gridNumCells[2], /*gridSize,*/ gridIndex[2], gridCellParticleListStart;
    int gridIndexMin[2], gridIndexMax[2];
    int i, j, k, s, t, p1, p2, total;
    
    /* Work out grid dimensions and allocate. */
    gridNumCells[0] = (int)(sqrt(maxParticles) + 1);
    gridNumCells[1] = (int)(sqrt(maxParticles) + 1);
    gridCellSize[0] = (arena.max[0] - arena.min[0]) / gridNumCells[0];
    gridCellSize[1] = (arena.max[1] - arena.min[1]) / gridNumCells[1];
    //gridSize = gridNumCells[0] * gridNumCells[1];
    
    /* Assumption. */
    for (i = 0; i < maxParticles; i++)
    if (particle[i].radius * 2.0 > gridCellSize[0] ||
        particle[i].radius * 2.0 > gridCellSize[1])
    panic("collideParticlesUniformGrid: particle diameter > cellSize\n");
    
    /* Allocate arrays. */
    gridCellParticleCount = (int **)malloc(gridNumCells[0] * sizeof(int *));
    if (gridCellParticleCount == 0)
    panic("collideParticlesUniformGrid: malloc failed\n");
    gridCellParticleListEnd = (int **)malloc(gridNumCells[0] * sizeof(int *));
    if (gridCellParticleListEnd == 0)
    panic("collideParticlesUniformGrid: malloc failed\n");
    for (i = 0; i < gridNumCells[0]; i++) {
        gridCellParticleCount[i] = (int *)malloc(gridNumCells[1] * sizeof(int));
        if (gridCellParticleCount[i] == 0)
        panic("collideParticlesUniformGrid: malloc failed\n");
        gridCellParticleListEnd[i] = (int *)malloc(gridNumCells[1] * sizeof(int));
        if (gridCellParticleListEnd[i] == 0)
        panic("collideParticlesUniformGrid: malloc failed\n");
    }
    gridCellParticleList = (int *)malloc(maxParticles * sizeof(int));
    
    /* Initialise grid particle count. */
    for (i = 0; i < gridNumCells[0]; i++)
    for (j = 0; j < gridNumCells[1]; j++)
    gridCellParticleCount[i][j] = 0;
    
    /* Cell counts. */
    for (i = 0; i < maxParticles; i++) {
        calcGridIndex(particle[i], arena, gridCellSize, gridNumCells, gridIndex);
        gridCellParticleCount[gridIndex[0]][gridIndex[1]] += 1;
    }
    
    if (debug[debug_collideParticlesUniformGrid]) {
        printf("collideParticlesUniformGrid: gridCellParticleCount\n");
        for (i = 0; i < gridNumCells[0]; i++)
        for (j = 0; j < gridNumCells[1]; j++)
        printf("%d %d %d\n", i, j, gridCellParticleCount[i][j]);
    }
    
    /* Work out end of cell lists by accumulating cell counts. */
    for (i = 0; i < gridNumCells[0]; i++)
    for (j = 0; j < gridNumCells[1]; j++)
    gridCellParticleListEnd[i][j] = 0;
    
    total = 0;
    for (i = 0; i < gridNumCells[0]; i++)
    for (j = 0; j < gridNumCells[1]; j++) {
        total = total + gridCellParticleCount[i][j];
        gridCellParticleListEnd[i][j] = total - 1;
    }
    
    if (debug[debug_collideParticlesUniformGrid]) {
        printf("collideParticlesUniformGrid: gridCellParticleListEnd\n");
        for (i = 0; i < gridNumCells[0]; i++)
        for (j = 0; j < gridNumCells[1]; j++)
        printf("%d %d %d\n", i, j, gridCellParticleListEnd[i][j]);
    }
    
    /* Build particle lists. */
    for (i = 0; i < gridNumCells[0]; i++)
    for (j = 0; j < gridNumCells[1]; j++)
    gridCellParticleCount[i][j] = 0;
    
    for (i = 0; i < maxParticles; i++) {
        calcGridIndex(particle[i], arena, gridCellSize, gridNumCells, gridIndex);
        gridCellParticleList[gridCellParticleListEnd[gridIndex[0]][gridIndex[1]] -
                             gridCellParticleCount[gridIndex[0]][gridIndex[1]]] = i;
        gridCellParticleCount[gridIndex[0]][gridIndex[1]] += 1;
    }
    
    if (debug[debug_collideParticlesUniformGrid]) {
        printf("collideParticlesUniformGrid: gridCellParticleList\n");
        for (i = 0; i < gridNumCells[0]; i++) {
            for (j = 0; j < gridNumCells[1]; j++) {
                gridCellParticleListStart =
                gridCellParticleListEnd[i][j] - gridCellParticleCount[i][j] + 1;
                printf("particle list %d %d\n", i, j);
                for (k = gridCellParticleListStart;
                     k < gridCellParticleListEnd[i][j];
                     k++)
                printf("%d\n", gridCellParticleList[k]);
                printf("\n");
            }
        }
    }
    
    /* Collision detection. */
    for (i = 0; i < maxParticles; i++) {
        calcGridIndex(particle[i], arena, gridCellSize, gridNumCells, gridIndex);
        
        /* Grid index bounds for this particle. */
        gridIndexMin[0] = gridIndex[0] - 1;
        if (gridIndexMin[0] < 0)
        gridIndexMin[0] = 0;
        gridIndexMin[1] = gridIndex[1] - 1;
        if (gridIndexMin[1] < 0)
        gridIndexMin[1] = 0;
        gridIndexMax[0] = gridIndex[0] + 1;
        if (gridIndexMax[0] > gridNumCells[0] - 1)
        gridIndexMax[0] = gridNumCells[0] - 1;
        gridIndexMax[1] = gridIndex[1] + 1;
        if (gridIndexMax[1] > gridNumCells[1] - 1)
        gridIndexMax[1] = gridNumCells[1] - 1;
        
        p1 = i;
        
        for (s = gridIndexMin[0]; s <= gridIndexMax[0]; s++) {
            for (t = gridIndexMin[1]; t <= gridIndexMax[1]; t++) {
                gridCellParticleListStart =
                gridCellParticleListEnd[s][t] - gridCellParticleCount[s][t] + 1;
                for (j = gridCellParticleListStart;
                     j <= gridCellParticleListEnd[s][t];
                     j++) {
                    p2 = gridCellParticleList[j];
                    
                    /* Don't test particle against itself. */
                    if (p2 == p1)
                    continue;
                    
                    /* Only test pairs once. */
                    if (p2 < p1)
                    continue;
                    
                    if (debug[debug_collideParticlesUniformGrid])
                    printf("collideParticlesUniformGrid: testing %d %d\n", p1, p2);
                    
                    if (collisionDetectionParticles(particle[p1], particle[p2])) {
                        
                        if (debug[debug_collideParticlesUniformGrid])
                        printf("collision: %d %d\n", p1, p2);
                        
                        /* Take step back. Better approaches possible. */
                        eulerStepSingleParticle(particle[p1], -h);
                        eulerStepSingleParticle(particle[p2], -h);
                        
                        if (debug[debug_collideParticlesUniformGrid]) {
                            printf("velocities before: %f %f %f %f\n",
                                   particle[p1].velocity[0], particle[p1].velocity[1],
                                   particle[p2].velocity[0], particle[p2].velocity[1]);
                        }
                        
                        /* Collision */
                        if (dimension == 1)
                        collisionReactionParticles1D(particle[p1], particle[p2]);
                        else if (dimension == 2) {
                            if (reacCalc == basisChange)
                            collisionReactionParticles2DbasisChange(particle[p1], particle[p2]);
                            else if (reacCalc == projNormal)
                            collisionReactionParticles2DprojNormal(particle[p1], particle[p2]);
                            else
                            panic("collision reaction calculation not specified\n");
                        }
                        
                        if (debug[debug_collideParticlesUniformGrid])
                        printf("velocities after: %f %f %f %f\n",
                               particle[p1].velocity[0], particle[p1].velocity[1],
                               particle[p1].velocity[0], particle[p2].velocity[1]);
                        
                        /* Step forward. */
                        eulerStepSmoothAbsorption(particle[p1], h);
                        eulerStepSmoothAbsorption(particle[p2], h);
                        
                        /* Check walls. */
                        collideParticleWall(particle[p1], arena);
                        collideParticleWall(particle[p2], arena);
                    }
                }
            }
        }
    }
    
    /* Free arrays. */
    for (i = 0; i < gridNumCells[0]; i++) {
        free(gridCellParticleCount[i]);
        free(gridCellParticleListEnd[i]);
    }
    free(gridCellParticleCount);
    free(gridCellParticleListEnd);
    free(gridCellParticleList);
}

void updateParticles(void)
{
    static Real time = 0.0, h;
    
    if (!g.pause) {
        /* Calculate time increment. */
        h = g.elapsedTime - time;
        time = g.elapsedTime;
        if (debug[debug_time])
        printf("updateParticles: time %f %f\n", time, h);
    } else {
        h = 0.0;
    }
    
    /* Compute new positions of particles. */
    integrateMotionParticles(h);
    
    /* Collisions against walls. */
    collideParticlesWall();
    
    /* Collisions amongst particles. */
    if (CDmethod == bruteForce)
    collideParticlesBruteForce(h);
    else if (CDmethod == uniformGrid)
    collideParticlesUniformGrid(h);
    else
    panic("updateParticles: unknown collision detection method\n");
    if (g.gameLevel == 1) {
        if (g.numMaxRadius == 1) {
            if (particle[0].radius == g.maxRadius) {
                printf("Player wins!\nPress [SPACE] bar to Continue...\n");
                g.win = true;
                g.go = false;
            } else {
                printf("Player loses!\nPress [SPACE] bar to Restart the level...\n");
                g.lose = true;
                g.go = false;
            }
        } else
        g.numMaxRadius = 0;
    }
    
    if (debug[debug_sum_kinetic_energy]) {
        printf("K = %f\n", sumKineticEnergy());
    }
    if (debug[debug_sum_momentum]) {
        Real p[2];
        sumMomentum(p);
        printf("p = %f %f\n", p[0], p[1]);
    }
}

void displayParticles(void)
{
    /* to debug player, turn off the other particles for now */
    for (int i = 1; i < maxParticles; i++) {
        if (debug[debug_particle])
        printf ("displayParticles: x %f y %f\n",
                particle[i].position[0], particle[i].position[1]);
        glPushMatrix();
        // glColor3f (0.8, 0.8, 0.8);
        glTranslatef(particle[i].position[0], particle[i].position[1], 0.0);
        displayParticle(&particle[i], 1.0, 1.0, 1.0);
        glPopMatrix();
    }
}

void displayPlayer(void)
{
    glPushMatrix();
    glColor3f(0.8, 0.0, 0.0);
    glTranslatef(particle[0].position[0], particle[0].position[1], 0.0);
    displayParticle(&particle[0], 1.0, 1.0, 1.0);
    glPopMatrix();
}

void playerEjectionReaction(Particle &p)
{
    particle[0].mass -= p.mass;
    particle[0].velocity[0] -= p.mass / particle[0].mass * p.velocity[0];
    particle[0].velocity[1] -= p.mass / particle[0].mass * p.velocity[1];
    particle[0].radius = sqrt(particle[0].mass);
}

#if 1
void ejectingMotes()
{
    int w, h;
    SDL_GetWindowSize(mainWindow, &w, &h);
    
    // filter noise values
    if (g.x > 10.0 && g.y > 10.0) {
        g.x = arena.min[0] + (Real)(g.x / w * (arena.max[0] - arena.min[0]));
        g.y = arena.min[1] + (Real)(g.y / h * (arena.max[1] - arena.min[1]));
        
        // ejecting motes to move player mote
        Real v = 10.0;
        Real alpha;
        alpha = 2.0 * M_PI - atan2f(g.y, g.x);
        
        int i = 1;
        while (!particle[i].isAbort && i < maxParticles)
        i++;
        
        particle[i].position[0] = particle[0].position[0] + (particle[0].radius + 0.1) * cosf(alpha);
        particle[i].position[1] = particle[0].position[1] + (particle[0].radius + 0.1) * sinf(alpha);
        particle[i].velocity[0] = v * cosf(alpha);
        particle[i].velocity[1] = v * sinf(alpha);
        particle[i].mass = 0.0002;
        particle[i].radius = sqrt(particle[i].mass);
        particle[i].elasticity = 1.0;
        particle[i].quadric = gluNewQuadric();
        particle[i].slices = 10;
        particle[i].loops = 3;
        particle[i].isAbort = false;
        playerEjectionReaction(particle[i]);
        if (particle[0].radius < 0.001) {
            g.lose = true;
            g.go = false;
            printf("Player lose!\nPress [SPACE] bar to Restart the level...\n");
        }
    }
}
#endif

void display()
{
    static int frameNo = 0;
    GLenum err;
    
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    glColor3f (0.8, 0.8, 0.8);
    
    glPushMatrix();
    
    /* Display particle and arena. */
    glScalef(camera.scale, camera.scale, 1.0);
    glTranslatef(-particle[0].position[0], -particle[0].position[1], 0.0);
    
    glPushMatrix();
    
    displayArena();
    if (g.renMode == solid)
    glDisable(GL_LIGHTING);
    glUseProgram(shaderProgram);
    
    displayPlayer();
    displayParticles();
    
    glUseProgram(0);
    if (g.renMode == solid)
    glEnable(GL_LIGHTING);
    
    glPopMatrix();
    
    glTranslatef(particle[0].position[0], particle[0].position[1], 0.0);
    glScalef(1.0 / camera.scale, 1.0 / camera.scale, 1.0);
    
    frameNo++;
    
    glPopMatrix();
    
    SDL_GL_SwapWindow(mainWindow);
    
    // Check for errors.
    while ((err = glGetError()) != GL_NO_ERROR) {
        printf("display: %s\n",gluErrorString(err));
    }
    
}

void update()
{
    if (!g.go)
    return;
    
    g.elapsedTime = SDL_GetTicks() / (Real)milli
    - g.startTime - g.pauseInterval;
    
    if (g.eject)
    ejectingMotes();
    
    updateParticles();
    
    postRedisplay();
}

bool handleMouseMotion(SDL_MouseMotionEvent e)
{
    float dy;
    
    if (debug[debug_mouse]) {
        printf("motion: %d %d\n", e.x, e.y);
        printf("camera.scale: %f\n", camera.scale);
    }
    
    dy = (float)(e.y - camera.lastY);
    camera.lastY = e.y;
    g.x = e.x;
    g.y = e.y;
    
    switch (camera.control) {
        case inactive:
        break;
        case zoom:
        camera.scale += dy/100.0;
        break;
    }
    
    return false;
}

bool handleMouseButtonDown(SDL_MouseButtonEvent e)
{
    if (g.go) {
        camera.lastY = e.y;
        if (!g.pause) {
            g.x = e.x;
            g.y = e.y;
        }
    }
#if 1
    if (e.button == SDL_BUTTON_LEFT)
    g.eject = true;
#endif
    if (e.button == SDL_BUTTON_RIGHT)
    camera.control = zoom;
    
    return false;
}

bool handleMouseButtonUp(SDL_MouseButtonEvent e)
{
#if 1
    if (e.button == SDL_BUTTON_LEFT)
    g.eject = false;
#endif
    camera.control = inactive;
    
    return false;
}

bool handleKeyDown(SDL_Keysym *key)
{
    // Handle key press events here
    switch (key->sym)
    {
        case SDLK_ESCAPE:
        return true;
        case SDLK_w:
        changeRenderMode();
        postRedisplay();
        break;
        case SDLK_r:
        if (reacCalc == projNormal)
        reacCalc = basisChange;
        else if (reacCalc == basisChange)
        reacCalc = projNormal;
        break;
        case SDLK_SPACE:
        if (!g.go) {
            g.startTime = SDL_GetTicks() / (Real)milli;
            g.go = true;
            if (g.win) {
                if (g.gameLevel < MAX_GAME_LEVEL)
                g.gameLevel++;
                initStage(g.gameLevel);
                g.win = false;
            }
            if (g.lose) {
                g.lose = false;
                initStage(g.gameLevel);
            }
        } else {
            g.pause = !g.pause;
            if (g.pause)
            g.pauseTime = SDL_GetTicks() / (Real)milli;
            else {
                g.resumeTime = SDL_GetTicks() / (Real)milli;
                g.pauseInterval += g.resumeTime - g.pauseTime;
            }
        }
        break;
        default:
        break;
    }
    
    return false;
}

bool handleEvents()
{
    SDL_Event event;
    
    while (SDL_PollEvent(&event))
    {
        switch (event.type)
        {
            case SDL_MOUSEMOTION:
            return handleMouseMotion(event.motion);
            break;
            
            case SDL_MOUSEBUTTONDOWN:
            return handleMouseButtonDown(event.button);
            break;
            
            case SDL_MOUSEBUTTONUP:
            return handleMouseButtonUp(event.button);
            break;
            
            case SDL_KEYDOWN:
            return handleKeyDown(&event.key.keysym);
            break;
            
            case SDL_WINDOWEVENT:
            switch (event.window.event)
            {
                case SDL_WINDOWEVENT_RESIZED:
                if (event.window.windowID == SDL_GetWindowID(mainWindow))
                SDL_SetWindowSize(mainWindow, event.window.data1, event.window.data2);
                break;
                case SDL_WINDOWEVENT_CLOSE:
                return true;
                break;
                default:
                break;
            }
            break;
            
            default:
            break;
        }
    }
    
    return false;
}

void mainLoop()
{
    while (1) {
        if (handleEvents()) {
            break;
        }
        if(wantRedisplay) {
            display();
            wantRedisplay = 0;
        }
        update();
    }
}

void reshape (int w, int h)
{
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    
    glOrtho(-10.0, 10.0, -10.0, 10.0, -10.0, 10.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}


int main(int argc, char **argv)
{
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        fprintf(stderr, "Unable to init SDL: %s\n", SDL_GetError());
        return EXIT_SUCCESS;
    }
    
    SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 8);
    SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 8);
    SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 8);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    
    mainWindow = SDL_CreateWindow("Assignment 3 Osmos-like Game",
                                  500, 500, 600, 600, SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE);
    if (!mainWindow) {
        fprintf(stderr, "Failed to create a window: %s\n", SDL_GetError());
        return EXIT_SUCCESS;
    }
    
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 2);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 1);
    SDL_GLContext mainGLContext = SDL_GL_CreateContext(mainWindow);
    if (mainGLContext == 0) {
        fprintf(stderr, "Unable to get OpenGL context: %s\n", SDL_GetError());
        return EXIT_SUCCESS;
    }
    
    if (SDL_GL_MakeCurrent(mainWindow, mainGLContext) != 0) {
        fprintf(stderr, "Unable to make OpenGL context current: %s\n", SDL_GetError());
        return EXIT_SUCCESS;
    }
    
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 2);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 1);
    int w, h;
    SDL_GetWindowSize(mainWindow, &w, &h);
    reshape(w, h);
    
    init();
    
    // Main event and display loop goes here
    mainLoop();
    
    SDL_DestroyWindow(mainWindow);
    SDL_Quit();
    
    return EXIT_SUCCESS;
}
