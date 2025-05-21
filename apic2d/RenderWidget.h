#ifndef RENDER_WIDGET_H
#define RENDER_WIDGET_H


#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "imgui.h"
#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"

#include <chrono>
#include <vector>
#include "Shader.h"
#include <glm/glm.hpp>
#include "ParticalSystem.h"

class RenderWidget
{
public:
    RenderWidget();
    ~RenderWidget();

    int32_t Init();

    void Update();

    int32_t Destroy();

    bool ShouldClose();

    void ProcessInput();

    void PollEvents();

    void LoadVertexes(Fluid2d::ParticalSystem& ps);

    void LoadVertexes_new(std::vector<glm::vec2>& mPositions, std::vector<float>& mDensity);
    
    void DrawImGuiSidebar();

private:
    bool CreateRenderWindow(); // 原名称CreateWindow似乎会与原项目的其他库函数冲突

    float CalculateFPS();

    static void ResizeCallback(GLFWwindow* window, int width, int height);

private:

    GLFWwindow* mWindow = nullptr;
    int mWindowWidth = 1250;
    int mWindowHeight = 1000;

    Glb::Shader* mParticalShader = nullptr;
    Glb::Shader* mSdfShader = nullptr;
    Glb::Shader* mMilkShader = nullptr;

    GLuint mVaoParticals = 0;
    GLuint mPositionBuffer = 0;
    GLuint mDensityBuffer = 0;

    GLuint mFboSdf = 0;
    GLuint mTextureSdf = 0;
    GLuint mRboSdf = 0;

    size_t mParticalNum = 0;

    std::chrono::system_clock::time_point mUpdateTime;
    
    
};


#endif // !RENDER_WIDGET_H

