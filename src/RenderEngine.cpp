#include <CanvasPoint.h>
#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Loader.h>
#include <ModelTriangle.h>
#include <RayTriangleIntersection.h>
#include <TextureMap.h>

#include <chrono>
#include <cmath>
#include <map>
#include <random>
#include <string>
#include <vector>

#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>

enum Mode {
    RASTER,
    RAY,
    WIREFRAME
} MODE;

struct Camera {
    glm::mat4 transform;
    float near;
    float far;
    float fov;
};

struct Light {
    glm::mat4 transform;
    float radius;
    glm::vec3 colour;
};

// Standard normal distribution
std::default_random_engine gen;
std::normal_distribution<float> n_dist(0.0f, 1.0f); 

// Store all currently loaded texture maps
std::map<std::string, TextureMap> loaded_textures;

// Movement variables
bool rotate = false;


/*****************************
 * HELPERS
******************************/

// Switch rendering mode
void toggle_mode() {
    if (MODE == WIREFRAME) {
        MODE = RASTER;
    } else if (MODE == RASTER) {
        MODE = RAY;
    } else {
        MODE = WIREFRAME;
    }
}

// Get packed colour from colour vec3
uint32_t get_colour(glm::vec3 col) {
    return (uint32_t) (255 << 24) + ((int)col.x << 16) + ((int)col.y << 8) + (int)col.z;
}

// Get packed colour from colour vec3
glm::vec3 get_rgb(uint32_t col) {
    float b = 0xFF & col;
    float g = 0xFF & (col >> 8);
    float r = 0xFF & (col >> 16);
    return glm::vec3(r, g, b);
}

// Clamp
glm::vec3 clamp(glm::vec3 c, float max, float min) {
    if (c.x > max) c.x = max;
    if (c.y > max) c.y = max;
    if (c.z > max) c.z = max;
    if (c.x < min) c.x = min;
    if (c.y < min) c.y = min;
    if (c.z < min) c.z = min;
    return c;
}

// Apply a brightness to a colour
glm::vec3 intensify(glm::vec3 c, float b) {
    b = b < 0 ? 0 : b;
    c *= b;
    return clamp(c, 255, 0);
}

// All the barycentric interpolation things
template<typename T>
T interpolate_bary(T a, T b, T c, glm::vec3 w) {
    return a * w.x + b * w.y + c * w.z;
}

// Create a perspective projection matrix, for a frustrum
glm::mat4 perspective_projection(float fov, float aspect_ratio, float near, float far) {
    glm::mat4 proj(0);
    float s = 1 / tan(fov * M_PI/360);
    proj[0][0] = s / aspect_ratio;
    proj[1][1] = s;
    proj[2][2] = (far + near) / (near - far);
    proj[2][3] = -1.0f;
    proj[3][2] = 2 * far * near / (near - far);
    return proj;
}

// Create an orthographic projection matrix
glm::mat4 orthographic_projection(float fov, float aspect_ratio, float near, float far) {
    glm::mat4 proj(0);
    float s = 1 / tan(fov * M_PI/360);
    proj[0][0] = s / aspect_ratio;
    proj[1][1] = s;
    proj[2][2] = (far + near) / (near - far);
    proj[3][3] = 1.0f;
    proj[3][2] = 2 * far * near / (near - far);
    return proj;
}

/*********************
OBJECT TRANSFORMATIONS
**********************/

void rotate_x(glm::mat4 &transform, float theta) {
    glm::mat3 x_rot(1, 0, 0, 
                    0, cos(theta), sin(theta),
                    0, -sin(theta), cos(theta));
    transform *= glm::mat4(x_rot);
}

void rotate_y(glm::mat4 &transform, float theta) {
    glm::mat3 y_rot(cos(theta), 0, -sin(theta), 
                    0, 1, 0,
                    sin(theta), 0, cos(theta));
    transform *= glm::mat4(y_rot);
}

void orbit_camera_y(Camera &c, float theta) {
    glm::mat3 y_rot(cos(theta), 0, -sin(theta), 
                    0, 1, 0,
                    sin(theta), 0, cos(theta));
    c.transform[3] = glm::mat4(y_rot) * c.transform[3];
}

// Force the camera to rotate itself to look at a given position
void camera_look_at(Camera &c, glm::vec4 pos) {
    glm::vec3 forward = glm::normalize(glm::vec3(c.transform[3] - pos));
    glm::vec3 right = glm::cross(glm::vec3(0, 1, 0), forward);
    glm::vec3 up = glm::cross(forward, right);
    c.transform[0] = glm::vec4(right, 0);
    c.transform[1] = glm::vec4(up, 0);
    c.transform[2] = glm::vec4(forward, 0);
}

// Simply rotates camera around the scene
void orbit(Camera &camera, float amount) {
    orbit_camera_y(camera, amount);
    camera_look_at(camera, glm::vec4(0, 0, 0, 1));
}


/******************
RAYTRACING HELPERS
******************/

// Find the closest intersection from a point with a direction 
RayTriangleIntersection get_closest_intersection(glm::vec3 start, glm::vec3 dir, std::vector<ModelTriangle> triangles) {
    RayTriangleIntersection closest(glm::vec3(0), 999999999999, ModelTriangle(), -1);
    for (int i = 0; i < triangles.size(); i++) {
        ModelTriangle triangle = triangles[i];
        // Calculate the triangle plane
        glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
        glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
        glm::vec3 relative_dir = start - triangle.vertices[0];
        glm::mat3 plane(-dir, e0, e1);

        glm::vec3 point = glm::inverse(plane) * relative_dir;
        if (point.x >= 0 && point.y >= 0.0 && point.y <= 1.0 && point.z >= 0.0 && point.z <= 1.0 && point.y + point.z <= 1.0) { 
            if (point.x < closest.distance) {
                closest = RayTriangleIntersection(point, point.x, triangle, i);
            }
        }
    }
    return closest;
}

// Check if object is in shadow
bool in_shadow(glm::vec3 start, glm::vec3 dir_to_light, std::vector<ModelTriangle> triangles, float dist_to_light) {
    for (int i = 0; i < triangles.size(); i++) {
        ModelTriangle triangle = triangles[i];
        // Calculate the triangle plane
        glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
        glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
        glm::vec3 relative_dir = start - triangle.vertices[0];
        glm::mat3 plane(-dir_to_light, e0, e1);

        glm::vec3 point = glm::inverse(plane) * relative_dir;
        if (point.x >= 0 && point.x <= dist_to_light && point.y >= 0.0 && point.y <= 1.0 && point.z >= 0.0 && point.z <= 1.0 && point.y + point.z <= 1.0) {
            return true;
        }
    }
    return false;
}


/*****************************
 * DRAWING
******************************/

// Draw a pixel to the screen taking into account its depth relative to the camera
void draw_pixel(DrawingWindow &window, std::vector<std::vector<float>> &depth_buffer, int x, int y, uint32_t col, float depth) {
    if (x >= window.width || x < 0 || y >= window.height || y < 0 || depth <= 0) return;
    if (depth_buffer[y][x] > depth || depth_buffer[y][x] == 0) {
        depth_buffer[y][x] = depth;
        window.setPixelColour(x, y, col);
    }
}

// Simple DDA line drawing algorithm with colour
void draw_line(DrawingWindow &window, CanvasPoint start, CanvasPoint end) {
    float dx = end.x - start.x;
    float dy = end.y - start.y;
    glm::vec3 dc = end.colour - start.colour;

    float points = std::max(abs(dx), abs(dy));
    float inv_points = 1 / points;
    float xi = dx * inv_points;
    float yi = dy * inv_points;
    glm::vec3 ci = dc * inv_points;

    float x = start.x;
    float y = start.y;
    glm::vec3 c = start.colour;

    for (int i = 0; i < points; i++) {
        window.setPixelColour(floor(x), floor(y), get_colour(c));
        x += xi;
        y += yi;
        c += ci;
    }
}

// Draw an unfilled triangle
void draw_stroked_triangle(DrawingWindow &window, CanvasTriangle triangle) {
    draw_line(window, triangle.v0(), triangle.v1());
    draw_line(window, triangle.v0(), triangle.v2());
    draw_line(window, triangle.v1(), triangle.v2());
}

// Use barycentric filling to rasterize a triangle
void draw_filled_triangle(DrawingWindow &window, std::vector<std::vector<float>> &depth_buffer, CanvasTriangle triangle, TextureMap texture = TextureMap()) {
    CanvasPoint v0 = triangle.v0();
    CanvasPoint v1 = triangle.v1();
    CanvasPoint v2 = triangle.v2();

    if (v0.depth < 0 && v1.depth < 0 && v2.depth < 0) return;

    // Calculate bounding box
    int minX = std::min({v0.x, v1.x, v2.x});
    int minY = std::min({v0.y, v1.y, v2.y});
    int maxX = std::max({v0.x, v1.x, v2.x});
    int maxY = std::max({v0.y, v1.y, v2.y});

    // Clip the values so they are within the screen
    minX = std::max(minX, 0);
    minY = std::max(minY, 0);
    maxX = std::min(maxX, window.width - 1);
    maxY = std::min(maxY, window.height - 1);

    // Check if each pixel in the bounding box is within the triangle or not
    float inv_area = 1.0f / ((v2.x - v0.x) * (v1.y - v0.y) - (v2.y - v0.y) * (v1.x - v0.x));
    for (int y = minY; y <= maxY; y++) {
        for (int x = minX; x <= maxX; x++) {
            glm::vec3 pixel_colour = glm::vec3(0, 0, 0);
            float pixel_depth = 0.0f;

            // Compute barycentric coords of pixel centre
            glm::vec2 p(x + 0.5f, y + 0.5f);
            float u = (p.x - v1.x) * (v2.y - v1.y) - (p.y - v1.y) * (v2.x - v1.x);
            float v = (p.x - v2.x) * (v0.y - v2.y) - (p.y - v2.y) * (v0.x - v2.x);
            float w = (p.x - v0.x) * (v1.y - v0.y) - (p.y - v0.y) * (v1.x - v0.x);

            bool inside = true;
            inside &= u == 0 ? (v2.y - v1.y == 0 && v2.x - v1.x > 0) || (v2.y - v1.y > 0) : u > 0;
            inside &= v == 0 ? (v0.y - v2.y == 0 && v0.x - v2.x > 0) || (v0.y - v2.y > 0) : v > 0;
            inside &= w == 0 ? (v1.y - v0.y == 0 && v1.x - v0.x > 0) || (v1.y - v0.y > 0) : w > 0;

            // Test if pixel is in the triangle
            if (inside) {
                glm::vec3 bary = glm::vec3(u, v, w) * inv_area; 
                // Calculate the depth at this pixel
                pixel_depth = interpolate_bary(v0.depth, v1.depth, v2.depth, bary);
                
                // Check if the triangle is textured
                if (v0.texturePoint == glm::vec2(-1,-1)) {
                    pixel_colour = interpolate_bary(v0.colour, v1.colour, v2.colour, bary);
                } else {
                    glm::vec2 textureCoords = interpolate_bary(v0.texturePoint, v1.texturePoint, v2.texturePoint, bary) / pixel_depth;
                    uint32_t colour = texture.pixels[round(textureCoords.x) + round(textureCoords.y) * texture.width];
                    pixel_colour = get_rgb(colour);
                }

                draw_pixel(window, depth_buffer, x, y, get_colour(pixel_colour), pixel_depth);
            }
        }
    }
}

// Use raytracing to render a more realistic scene at the cost of speed
void draw_raytraced(DrawingWindow &window, Camera camera, Light light, std::vector<ModelTriangle> model) {
    // Calculate pixel scale values
    float py_scale = tan((camera.fov * M_PI / (float) 360));
    float px_scale = py_scale * window.width / (float) window.height; 
    
    // Prepare values for anti_aliasing
    float aa_size = 2;
    float aa_samples = aa_size * aa_size;
    float aa_inc = 1 / aa_size; 
    float aa_offset = aa_inc / 2.0f; 
    for (int y = 0; y < window.height; y++) {
        for (int x = 0; x < window.width; x++) {
            
            glm::vec3 pixel_colour;
            float yn = y + aa_offset;
            for (int l = 0; l < aa_size; l++) {
                float xn = x + aa_offset;
                for (int k = 0; k < aa_size; k++) {
                    // Find the pixel coordinates on the image plane
                    float py = (1 - (2 * yn - 1) / window.height) * py_scale;
                    float px = ((2 * xn + 1) / window.width - 1) * px_scale;


                    // Find the object that intersects the pixel ray
                    glm::vec3 dir_to_pixel = glm::vec3(glm::normalize(camera.transform * glm::vec4(px, py, -1, 0)));
                    RayTriangleIntersection intersection = get_closest_intersection(glm::vec3(camera.transform[3]), dir_to_pixel, model);
                    if (intersection.triangle_index != -1) {
                        // Compute intersection position
                        glm::vec3 pos = glm::vec3(camera.transform[3]) + intersection.point.x * dir_to_pixel;
                        glm::vec3 bary(1 - intersection.point.y - intersection.point.z, intersection.point.y, intersection.point.z);

                        // Compute the surface normal of the point
                        glm::vec3 point_normal;
                        if (intersection.triangle.has_vns) {
                            point_normal = glm::normalize(interpolate_bary(intersection.triangle.vertex_normals[0], intersection.triangle.vertex_normals[1], intersection.triangle.vertex_normals[2], bary));
                        } else {
                            point_normal = intersection.triangle.surface_normal;
                        }
                        // Prevent acne
                        glm::vec3 acne_pos = pos + 0.0001f * point_normal;

                        // Check if the material is reflective
                        Material mat = intersection.triangle.material;
                        if (mat.model > 2) {
                            glm::vec3 reflected_ray = dir_to_pixel - point_normal * 2.0f * glm::dot(dir_to_pixel, point_normal);
                            RayTriangleIntersection reflected = get_closest_intersection(pos + 0.0001f * point_normal, glm::normalize(reflected_ray), model);
                            mat = reflected.triangle.material;
                        }

                        // Lighting with soft shadows
                        int shadow_samples = 1;
                        glm::vec3 colour;
                        for (int i = 0; i < shadow_samples; i++) {
                            // Get random point on light sphere
                            glm::vec3 sphere_point(n_dist(gen), n_dist(gen), n_dist(gen));
                            sphere_point = glm::normalize(sphere_point) * light.radius;

                            // Compute lighting effects
                            glm::vec3 vec_to_light = sphere_point + glm::vec3(light.transform[3]) - pos;
                            glm::vec3 dir_to_light = glm::normalize(vec_to_light);

                            // Calculate light if we are not in shadow
                            if (!in_shadow(acne_pos, dir_to_light, model, glm::length(vec_to_light))) {
                                // Proximity
                                float diffuse_intensity = 1 / (float) (4 * M_PI * glm::dot(vec_to_light, vec_to_light));
                                // Angle of incidence
                                diffuse_intensity += glm::dot(point_normal, dir_to_light);

                                // Specular
                                glm::vec3 reflected_ray = -dir_to_light + point_normal * 2.0f * glm::dot(dir_to_light, point_normal);
                                float specular_intensity = pow(std::max({glm::dot(glm::normalize(reflected_ray), -dir_to_pixel), 0.0f}), mat.shininess);

                                colour += intensify(mat.diffuse_colour, diffuse_intensity) + intensify(mat.specular_colour, specular_intensity);
                            }
                        }
                        
                        // Ambient
                        float ambient_intensity = 0.2f;

                        // Calculate final colour of pixel
                        colour /= (float) shadow_samples;
                        colour += intensify(mat.ambient_colour, ambient_intensity);
                    
                        pixel_colour += clamp(colour, 255, 0);
                    }
                    xn += aa_inc;
                }
                yn += aa_inc;
            }
            window.setPixelColour(x, y, get_colour(pixel_colour / aa_samples));
        }
    }
}

// Use rasterizing to render model, quicker approach that sacrifces on realism
void draw_raster(DrawingWindow &window, Camera camera, Light light, std::vector<ModelTriangle> model) {
    glm::mat4 proj = perspective_projection(camera.fov, window.width / (float) window.height, camera.near, camera.far);
    glm::mat4 ndc_to_canvas(1);
    ndc_to_canvas[0][0] = window.width / 2;
    ndc_to_canvas[1][1] = window.height / 2;
    glm::mat4 project_to_screen = ndc_to_canvas * proj * glm::inverse(camera.transform);

    std::vector<std::vector<float>> depth_buffer(window.height, std::vector<float>(window.width, 0));

    CanvasTriangle img_triangle;
    for (int j = 0; j < model.size(); j++) {
        ModelTriangle t = model[j];
    
        // Prevent backs of triangles from not being rendered
        if (glm::dot((t.vertices[0] - glm::vec3(camera.transform[3])), t.surface_normal) >= 0) {
            std::swap(t.vertices[2], t.vertices[1]);
        }
        // Project each vertex
        for (int i = 0; i < 3; i++) {
            glm::vec4 point = project_to_screen * glm::vec4(t.vertices[i], 1);
            img_triangle.vertices[i] = CanvasPoint(floor(point.x / point.w + window.width/2), floor(window.height/2 - point.y / point.w));
            img_triangle.vertices[i].depth = point.z / point.w;
            img_triangle.vertices[i].colour = t.material.diffuse_colour;
        }

        // Handle textures
        TextureMap tex;
        if (t.material.texture.compare("")) {
            if (loaded_textures.find(t.material.texture) == loaded_textures.end()) {
                tex = TextureMap("models/" + t.material.texture);
                loaded_textures[t.material.texture] = tex;
            } else {
                tex = loaded_textures[t.material.texture];
            }

            for (int i = 0; i < 3; i++) { 
                t.texturePoints[i].x *= tex.width;
                t.texturePoints[i].y *= tex.height;
                img_triangle.vertices[i].texturePoint = t.texturePoints[i] * img_triangle.vertices[i].depth;
            }
        }
        draw_filled_triangle(window, depth_buffer, img_triangle, tex);
    }
}

// Simply render the outlines of each object in the model
void draw_wireframe(DrawingWindow &window, Camera camera, std::vector<ModelTriangle> model) {
    glm::mat4 proj = perspective_projection(camera.fov, window.width / (float) window.height, camera.near, camera.far);
    glm::mat4 ndc_to_canvas(1);
    ndc_to_canvas[0][0] = window.width / 2;
    ndc_to_canvas[1][1] = window.height / 2;
    glm::mat4 project_to_screen = ndc_to_canvas * proj * glm::inverse(camera.transform);

    CanvasTriangle img_triangle;
    for (ModelTriangle t : model) {
        
        for (int i = 0; i < 3; i++) {
            glm::vec4 point = project_to_screen * glm::vec4(t.vertices[i], 1);
            img_triangle.vertices[i] = CanvasPoint(floor(point.x / point.w + window.width/2), floor(window.height/2 - point.y / point.w));
            img_triangle.vertices[i].depth = point.z / point.w;
            img_triangle.vertices[i].colour = glm::vec3(255,255,255);
        }
        draw_stroked_triangle(window, img_triangle);
    }
}

// Renders the cornell box scene from its obj file
void draw_model(DrawingWindow &window, Camera camera, Light light, std::vector<ModelTriangle> model) {
    // Set up timer
    auto t1 = std::chrono::high_resolution_clock::now();

    // Clear screen 
    window.clearPixels();

    std::string render_type;
    if (MODE == RAY) {
        render_type = "RAYTRACING";
        draw_raytraced(window, camera, light, model);
    } else if (MODE == RASTER) {
        render_type = "RASTERISING";
        draw_raster(window, camera, light, model);
    } else {
        render_type = "WIREFRAME";
        draw_wireframe(window, camera, model);
    }

    // Calculate time taken to render the model with the desired mode
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << render_type << " - Frame took " << duration << "ms" << std::endl;
}


/******************
 * SCENE CONTROL
******************/

void handleEvent(SDL_Event event, DrawingWindow &window, Camera &camera, Light &light, std::vector<ModelTriangle> model) {
    if (event.type == SDL_KEYDOWN) {
        if (event.key.keysym.sym == SDLK_LEFT) {
            if (rotate) {
                rotate_y(camera.transform, 0.1f);
            } else {
                camera.transform[3].x -= 0.1f;
            }
        } else if (event.key.keysym.sym == SDLK_RIGHT) {
            if (rotate) {
                rotate_y(camera.transform, -0.1f);
            } else {
                camera.transform[3].x += 0.1f;
            }
        } else if (event.key.keysym.sym == SDLK_UP) {
            if (rotate) {
                rotate_x(camera.transform, 0.1f);
            } else {
                camera.transform[3].y += 0.1f; 
            }
        } else if (event.key.keysym.sym == SDLK_DOWN) {
            if (rotate) {
                rotate_x(camera.transform, -0.1f);
            } else {
                camera.transform[3].y -= 0.1f; 
            }
        } else if (event.key.keysym.sym == SDLK_q) {
            camera.fov += 1;
        } else if (event.key.keysym.sym == SDLK_e) {
            camera.fov -= 1;
        } else if (event.key.keysym.sym == SDLK_TAB) {
            toggle_mode();
        } else if (event.key.keysym.sym == SDLK_r) {
            rotate = !rotate;
        }
        draw_model(window, camera, light, model);
    } else if (event.type == SDL_MOUSEBUTTONDOWN) window.savePPM("output.ppm");
}

int main(int argc, char *argv[]) {
    // Set up display window
    int width = 640;
    int height = 480;
    DrawingWindow window = DrawingWindow(width, height, false);
    SDL_Event event;

    // Setup camera, positive z out of screen
    glm::vec4 camera_pos(0, 0, 2, 1);
    glm::mat4 camera_transform(1);
    camera_transform[3] = camera_pos;
    Camera camera = {camera_transform, 0.1, 1000, 45};
    
    // Load in our model and choose the initial rendering mode
    std::vector<ModelTriangle> model = load_obj("textured-cornell-box", 0.17);
    MODE = WIREFRAME;

    // Create light source
    glm::vec4 light_pos(0, 0.35f, 0, 1);
    glm::mat4 light_transform(1);
    light_transform[3] = light_pos;
    Light light = {light_transform, 0.05f, glm::vec3(255,255,255)};

    draw_model(window, camera, light, model);
    while (true) {
        // We MUST poll for events - otherwise the window will freeze !
        if (window.pollForInputEvents(event)) handleEvent(event, window, camera, light, model);
        // Need to render the frame at the end, or nothing actually gets shown on the screen !
        window.renderFrame();
    }
}
