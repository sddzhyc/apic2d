#include "RenderWidget.h"

#include <iostream>
#include <fstream>


RenderWidget::RenderWidget() {

}

RenderWidget::~RenderWidget() {
    Destroy();
}

int32_t RenderWidget::Init() {
    if (!CreateWindow()) {
        return -1;
    }
    
    // ��ʼ��shader
    std::string particalVertShaderPath = "../code/Fluid2d/Shaders/DrawParticals.vert";
    std::string particalFragShaderPath = "../code/Fluid2d/Shaders/DrawParticals.frag";
    mParticalShader = new Glb::Shader();
    mParticalShader->BuildFromFile(particalVertShaderPath, particalFragShaderPath);

    std::string ballVertShaderPath = "../code/Fluid2d/Shaders/DrawSdf.vert";
    std::string ballGeomShaderPath = "../code/Fluid2d/Shaders/DrawSdf.geom";
    std::string ballFragShaderPath = "../code/Fluid2d/Shaders/DrawSdf.frag";
    mSdfShader = new Glb::Shader();
    mSdfShader->BuildFromFile(ballVertShaderPath, ballFragShaderPath, ballGeomShaderPath);

    std::string milkVertShaderPath = "../code/Fluid2d/Shaders/DrawMilk.vert";
    std::string milkFragShaderPath = "../code/Fluid2d/Shaders/DrawMilk.frag";
    mMilkShader = new Glb::Shader();
    mMilkShader->BuildFromFile(milkVertShaderPath, milkFragShaderPath);
    glUniform1i(glGetUniformLocation(mMilkShader->GetId(), "textureSdf"), 0);
    
    // ����vao
    glGenVertexArrays(1, &mVaoParticals);
    glGenBuffers(1, &mPositionBuffer);
    glGenBuffers(1, &mDensityBuffer);

    // ����SDF��frame buffer
    glGenFramebuffers(1, &mFboSdf);
    glBindFramebuffer(GL_FRAMEBUFFER, mFboSdf);
    // ��������
    glGenTextures(1, &mTextureSdf);
    glBindTexture(GL_TEXTURE_2D, mTextureSdf);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1000, 1000, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
    glBindTexture(GL_TEXTURE_2D, 0);
    // �󶨵�FBO
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, mTextureSdf, 0);
    // ����RBO
    glGenRenderbuffers(1, &mRboSdf);
    glBindRenderbuffer(GL_RENDERBUFFER, mRboSdf);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, 1000, 1000);
    glBindRenderbuffer(GL_RENDERBUFFER, 0);
    // �󶨵�FBO
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, mRboSdf);
    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
        std::cout << "ERROR: SDF Framebuffer is not complete!" << std::endl;
    }
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    // �ӿڴ�С
    glViewport(0, 0, mWindowWidth, mWindowHeight);
}

void RenderWidget::Update() {
    //// ������
    //glBindFramebuffer(GL_FRAMEBUFFER, 0);
    //glClearColor(0.3f, 0.3f, 0.3f, 1.0f);
    //glEnable(GL_DEPTH_TEST);
    //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    //glBindVertexArray(mVaoParticals);
    //mParticalShader->Use();
    //glEnable(GL_PROGRAM_POINT_SIZE);
    //glDrawArrays(GL_POINTS, 0, mParticalNum);

    // ��������
    glBindFramebuffer(GL_FRAMEBUFFER, mFboSdf);
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glBindVertexArray(mVaoParticals);
    mSdfShader->Use();
    glEnable(GL_PROGRAM_POINT_SIZE);
    glDrawArrays(GL_POINTS, 0, mParticalNum);

    // ��ţ��
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glClearColor(1.0f, 0.0f, 0.0f, 1.0f);
    glDisable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT);

    mMilkShader->Use();
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, mTextureSdf);
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

    // ��ʾFPS
    char title[128] = "";
    sprintf_s(title, "Fluid Simulation FPS=%.2f", CalculateFPS());
    glfwSetWindowTitle(mWindow, title);

    glfwSwapBuffers(mWindow);
}

int32_t RenderWidget::Destroy() {

    glDeleteVertexArrays(1, &mVaoParticals);
    glDeleteBuffers(1, &mPositionBuffer);
    glDeleteBuffers(1, &mDensityBuffer);
    delete mParticalShader;
    delete mSdfShader;
    delete mMilkShader;
    glfwTerminate();
    return 0;
}

bool RenderWidget::ShouldClose() {
    return glfwWindowShouldClose(mWindow);
}

void RenderWidget::ProcessInput() {
    if (glfwGetKey(mWindow, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
        glfwSetWindowShouldClose(mWindow, true);
    }

    return;
}

void RenderWidget::PollEvents() {
    glfwPollEvents();
}

void RenderWidget::LoadVertexes(Fluid2d::ParticalSystem& ps) {
    glBindVertexArray(mVaoParticals);
    glBindBuffer(GL_ARRAY_BUFFER, mPositionBuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec2) * ps.mPositions.size(), ps.mPositions.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, mDensityBuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * ps.mDensity.size(), ps.mDensity.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof(float), (void*)0);
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);

    mParticalNum = ps.mPositions.size();
}

bool RenderWidget::CreateWindow() {
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    

    // ��������
    mWindow = glfwCreateWindow(mWindowWidth, mWindowHeight, "Fluid Simulation", NULL, NULL);
    if (mWindow == nullptr)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return false;
    }
    glfwSetWindowPos(mWindow, 100, 100);
    glfwMakeContextCurrent(mWindow);

    // ע��ص�����
    glfwSetWindowUserPointer(mWindow, this);
    glfwSetFramebufferSizeCallback(mWindow, ResizeCallback);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return false;
    }

    return true;
}

float RenderWidget::CalculateFPS() {
    auto nowTime = std::chrono::system_clock::now();
    auto deltaTime = nowTime - mUpdateTime;
    mUpdateTime = nowTime;
    auto durMS = std::chrono::duration_cast<std::chrono::milliseconds>(deltaTime).count();
    float fps = 1000.0f / durMS;
    return fps;
}

void RenderWidget::ResizeCallback(GLFWwindow* window, int width, int height) {
    // �ҵ�thisָ��
    auto thisPtr = reinterpret_cast<RenderWidget*>(glfwGetWindowUserPointer(window));
    glViewport(0, 0, width, height);
}
