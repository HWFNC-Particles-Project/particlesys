#include <particle_array.h>

#ifdef NO_VIS
// dummy implementation
int particle_vis_init() {
    return 0;
}
void particle_vis_draw(particle_array *particles) {
    (void) particles;
}
void particle_vis_deinit() {
}
#else
#include <GL/glew.h>
#include <SDL2/SDL.h>
#include <stdio.h>
#include <stdlib.h>

static int check_shader_compile_status(GLuint obj) {
    GLint status;
    glGetShaderiv(obj, GL_COMPILE_STATUS, &status);
    if(status == GL_FALSE) {
        GLint length;
        glGetShaderiv(obj, GL_INFO_LOG_LENGTH, &length);
        char *log = malloc(length);
        glGetShaderInfoLog(obj, length, &length, &log[0]);
        fprintf(stderr, "%s", log);
        free(log);
        return 0;
    }
    return 1;
}

static int check_program_link_status(GLuint obj) {
    GLint status;
    glGetProgramiv(obj, GL_LINK_STATUS, &status);
    if(status == GL_FALSE) {
        GLint length;
        glGetProgramiv(obj, GL_INFO_LOG_LENGTH, &length);
        char *log = malloc(length);
        glGetProgramInfoLog(obj, length, &length, &log[0]);
        fprintf(stderr, "%s", log);
        free(log);
        return 0;
    }
    return 1;
}

typedef struct attribbind_t {
    const char *name;
    int location;
} attribbind;

static GLuint shader_build(const char *vertex_source, const char *fragment_source, const attribbind *attribs) {
    GLuint shader_program, vertex_shader, fragment_shader;

    vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex_shader, 1, &vertex_source, 0);
    glCompileShader(vertex_shader);
    if(!check_shader_compile_status(vertex_shader)) {
        return 0;
    }

    fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment_shader, 1, &fragment_source, 0);
    glCompileShader(fragment_shader);
    if(!check_shader_compile_status(fragment_shader)) {
        return 0;
    }

    shader_program = glCreateProgram();
    glAttachShader(shader_program, vertex_shader);
    glAttachShader(shader_program, fragment_shader);
    glDeleteShader(vertex_shader);
    glDeleteShader(fragment_shader);

    while(attribs->name != NULL) {
        glBindAttribLocation(shader_program, attribs->location, attribs->name);
        ++attribs;
    }

    glLinkProgram(shader_program);
    if(!check_program_link_status(shader_program)) {
        return 0;
    }
    return shader_program;
}
static void mul3(float *A, float *B, float *C) {
    for(int i = 0;i<4;++i) {
        for(int j = 0;j<4;++j) {
            float dot = 0;
            for(int k = 0;k<4;++k) {
                dot += B[i+4*k]*C[k+4*j];
            }
            A[i+4*j] = dot;
        }
    }
}

static void mul2(float *A, float *C) {
    float B[16];
    for(int i = 0;i<16;++i) {
        B[i] = A[i];
    }
    mul3(A,B,C);
}


void transform_identity(float *transform) {
    transform[0] = 1.0f; transform[4] = 0.0f; transform[8]  = 0.0f; transform[12] = 0.0f;
    transform[1] = 0.0f; transform[5] = 1.0f; transform[9]  = 0.0f; transform[13] = 0.0f;
    transform[2] = 0.0f; transform[6] = 0.0f; transform[10] = 1.0f; transform[14] = 0.0f;
    transform[3] = 0.0f; transform[7] = 0.0f; transform[11] = 0.0f; transform[15] = 1.0f;
}

void transform_scale(float *transform0, float sx, float sy, float sz) {
    float transform[16] = {
        sx,  0,  0,  0,
         0, sy,  0,  0,
         0,  0,  sz,  0,
         0,  0,  0,  1,
    };
    mul2(transform0, transform);
}

void transform_translate(float *transform0, float dx, float dy) {
    float transform[16] = {
         1,  0,  0,  0,
         0,  1,  0,  0,
         0,  0,  1,  0,
        dx, dy,  0,  1,
    };
    mul2(transform0, transform);
}

void transform_rotate_x(float *transform0, float angle) {
    float c = cos(angle);
    float s = sin(angle);
    float transform[16] = {
         1,  0,  0,  0,
         0,  c, -s,  0,
         0,  s,  c,  0,
         0,  0,  0,  1,
    };
    mul2(transform0, transform);
}

void transform_rotate_y(float *transform0, float angle) {
    float c = cos(angle);
    float s = sin(angle);
    float transform[16] = {
         c,  0, -s,  0,
         0,  1,  0,  0,
         s,  0,  c,  0,
         0,  0,  0,  1,
    };
    mul2(transform0, transform);
}

void transform_rotate_z(float *transform0, float angle) {
    float c = cos(angle);
    float s = sin(angle);
    float transform[16] = {
         c, -s,  0,  0,
         s,  c,  0,  0,
         0,  0,  1,  0,
         0,  0,  0,  1,
    };
    mul2(transform0, transform);
}

struct window_t {
    SDL_Window *sdlwindow;
    SDL_GLContext glcontext;
    GLuint shader_program;
    GLuint vbo;
}* window = NULL;

static const char *vertex_source =
    "#version 120\n"
    "uniform mat4 transform;\n"
    "attribute vec4 position;\n"
    "void main() {\n"
    "   gl_Position = transform*position;\n"
    "}\n";

static const char *fragment_source =
    "#version 120\n"
    "void main() {\n"
    "   gl_FragColor = vec4(1,1,1,1);\n"
    "}\n";

static const attribbind attribs[] = {
    {"position", 0},
    {NULL, 0}
};

static const int width = 640;
static const int height = 480;
static const size_t chunk_size = 1024;

int particle_vis_init() {
    window = malloc(sizeof(struct window_t));

    if(SDL_Init(SDL_INIT_VIDEO | SDL_INIT_EVENTS) < 0) {
        fprintf(stderr, "failed to init SDL");
        return 1;
    }

    if((window->sdlwindow = SDL_CreateWindow("Particle Visualization", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, width, height, SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN)) == 0) {
        fprintf(stderr, "failed to open window");
        SDL_Quit();
        free(window);
        window = NULL;
        return 1;
    }

    window->glcontext = SDL_GL_CreateContext(window->sdlwindow);

    if(glewInit()) {
        fprintf(stderr, "failed to init GLEW");
        SDL_GL_DeleteContext(window->glcontext);
        SDL_DestroyWindow(window->sdlwindow);
        SDL_Quit();
        free(window);
        window = NULL;
        return 1;
    }


    glGenBuffers(1, &window->vbo);
    glBindBuffer(GL_ARRAY_BUFFER, window->vbo);
    glBufferData(GL_ARRAY_BUFFER, chunk_size*4*sizeof(float), NULL, GL_STREAM_DRAW);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 4*sizeof(float), (char*)0 + 0*sizeof(GLfloat));

    window->shader_program = shader_build(vertex_source, fragment_source, attribs);
    glUseProgram(window->shader_program);
    return 0;
}

void particle_vis_draw(particle_array *particles) {
    (void) particles;
    SDL_Event event;
    while(SDL_PollEvent(&event)); // swallow events

    int mousex, mousey;
    SDL_GetMouseState(&mousex, &mousey);
    float alpha = 0.01*mousex;
    float beta = 0.01*mousey;

    glClearColor(0,0,0,1);
    glClear(GL_COLOR_BUFFER_BIT);
    float chunk[chunk_size*4];
    size_t size = particle_array_size(particles);

    float aspect = width/(float)height;
    float transform[16];
    transform_identity(transform);
    transform_scale(transform, 0.5f, 0.5f*aspect, 0.5f);
    transform_rotate_x(transform, -beta);
    transform_rotate_y(transform, alpha);

    glUniformMatrix4fv(glGetUniformLocation(window->shader_program, "transform"), 1, GL_FALSE, transform);

    for(size_t j = 0;j<size;j+=chunk_size) {
        size_t local_size = size-j;
        local_size = local_size<chunk_size?local_size:chunk_size;
        for(size_t i = 0;i<local_size;++i) {
            particle p = particle_array_get(particles, j+i);
            chunk[4*i+0] = p.position[0];
            chunk[4*i+1] = p.position[1];
            chunk[4*i+2] = p.position[2];
        }
        glBufferSubData(GL_ARRAY_BUFFER, 0, local_size*4*sizeof(float), chunk);
        glDrawArrays(GL_POINTS, 0, local_size);
    }

    SDL_GL_SwapWindow(window->sdlwindow);
}

void particle_vis_deinit() {
    if(window != NULL) {
        glDeleteProgram(window->shader_program);
        glDeleteBuffers(1, &window->vbo);
        SDL_GL_DeleteContext(window->glcontext);
        SDL_DestroyWindow(window->sdlwindow);
        SDL_Quit();
        free(window);
        window = NULL;
    }
}

#endif
