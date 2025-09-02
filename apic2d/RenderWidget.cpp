#include "RenderWidget.h"

#include <iostream>
#include <fstream>


RenderWidget::RenderWidget() {

}

RenderWidget::~RenderWidget() {
    Destroy();
}

int32_t RenderWidget::Init() {
    if (!CreateRenderWindow()) {
        return -1;
    }
    
    // 初始化shader
    //注意: APIC2D的项目工作目录是：D:\Documents\GitHub\apic2d\build\apic2d
    std::string particalVertShaderPath = "../../apic2d/Shaders/DrawParticals.vert";
    std::string particalFragShaderPath = "../../apic2d/Shaders/DrawParticals.frag";
    mParticalShader = new Glb::Shader();
    mParticalShader->BuildFromFile(particalVertShaderPath, particalFragShaderPath);

    std::string ballVertShaderPath = "../../apic2d/Shaders/DrawSdf.vert";
    std::string ballGeomShaderPath = "../../apic2d/Shaders/DrawSdf.geom";
    std::string ballFragShaderPath = "../../apic2d/Shaders/DrawSdf.frag";
    mSdfShader = new Glb::Shader();
    mSdfShader->BuildFromFile(ballVertShaderPath, ballFragShaderPath, ballGeomShaderPath);

    std::string milkVertShaderPath = "../../apic2d/Shaders/DrawMilk.vert";
    std::string milkFragShaderPath = "../../apic2d/Shaders/DrawMilk.frag";
    mMilkShader = new Glb::Shader();
    mMilkShader->BuildFromFile(milkVertShaderPath, milkFragShaderPath);
    glUniform1i(glGetUniformLocation(mMilkShader->GetId(), "textureSdf"), 0);
    
    // 生成vao
    glGenVertexArrays(1, &mVaoParticals);
    glGenBuffers(1, &mPositionBuffer);
    glGenBuffers(1, &mDensityBuffer);

    // 绘制SDF的frame buffer
    glGenFramebuffers(1, &mFboSdf);
    glBindFramebuffer(GL_FRAMEBUFFER, mFboSdf);
    // 生成纹理
    glGenTextures(1, &mTextureSdf);
    glBindTexture(GL_TEXTURE_2D, mTextureSdf);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1000, 1000, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
    glBindTexture(GL_TEXTURE_2D, 0);
    // 绑定到FBO
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, mTextureSdf, 0);
    // 生成RBO
    glGenRenderbuffers(1, &mRboSdf);
    glBindRenderbuffer(GL_RENDERBUFFER, mRboSdf);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, 1000, 1000);
    glBindRenderbuffer(GL_RENDERBUFFER, 0);
    // 绑定到FBO
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, mRboSdf);
    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
        std::cout << "ERROR: SDF Framebuffer is not complete!" << std::endl;
    }
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    // 视口大小
    glViewport(0, 0, mWindowWidth, mWindowHeight);
}

void RenderWidget::Update() {
  // 获取窗口尺寸
  int windowWidth, windowHeight;
  glfwGetWindowSize(mWindow, &windowWidth, &windowHeight);
  // 计算主内容区域的视口尺寸（假设侧边栏占窗口宽度的 20%）
  int mainViewportWidth = static_cast<int>(windowWidth * 0.8f);
  // 绑定默认帧缓冲并设置主内容视口
  glBindFramebuffer(GL_FRAMEBUFFER, 0);
  glViewport(0, 0, mainViewportWidth, windowHeight);  // 关键：限制主内容视口

  float rad0 = 40.0f;
  glm::vec2 minPos(-rad0, -rad0);
  glm::vec2 maxPos(rad0, rad0);
  glm::vec2 dataSize = maxPos - minPos;
      // 计算缩放比例（确保数据适应视口高度）
  float targetHeight = windowHeight;
  float targetWidth = mainViewportWidth;
  float scale = std::min(targetWidth / dataSize.x,  // 水平方向缩放比例
                         targetHeight / dataSize.y  // 垂直方向缩放比例
  );

  // 计算平移量（居中显示）
  glm::vec2 offset = -0.5f * (minPos + maxPos) * scale;          // 中心对齐
  offset += glm::vec2(0.5f * targetWidth, 0.5f * targetHeight);  // 视口中心

  // 设置视口和投影矩阵
  glViewport(0, 0, mainViewportWidth, windowHeight);
  glm::mat4 projection = glm::ortho(0.0f, static_cast<float>(mainViewportWidth),  // 左、右
                                    0.0f, static_cast<float>(windowHeight),       // 下、上
                                    -1.0f, 1.0f                                   // 近、远平面
  );

  // 应用缩放和平移（通过模型矩阵）
  glm::mat4 model = glm::mat4(1.0f);
  model = glm::translate(model, glm::vec3(offset, 0.0f));
  model = glm::scale(model, glm::vec3(scale, scale, 1.0f));
  // 将矩阵传递给着色器
  mParticalShader->Use();
  mParticalShader->SetMat4("projection", projection);
  mParticalShader->SetMat4("model", model);
  // 清屏并渲染
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glDrawArrays(GL_POINTS, 0, mParticalNum);
  // 清屏并渲染主内容（仅占用左侧 80% 区域）
  glClearColor(0.3f, 0.3f, 0.3f, 1.0f);
  glEnable(GL_DEPTH_TEST);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  // 画粒子
  glBindFramebuffer(GL_FRAMEBUFFER, 0);
  glClearColor(0.3f, 0.3f, 0.3f, 1.0f);
  glEnable(GL_DEPTH_TEST);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glBindVertexArray(mVaoParticals);
  mParticalShader->Use();
  glEnable(GL_PROGRAM_POINT_SIZE);
  glDrawArrays(GL_POINTS, 0, mParticalNum);

  // 画粒子球
  glBindFramebuffer(GL_FRAMEBUFFER, mFboSdf);
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  glEnable(GL_DEPTH_TEST);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glBindVertexArray(mVaoParticals);
  mSdfShader->Use();
  glEnable(GL_PROGRAM_POINT_SIZE);
  glDrawArrays(GL_POINTS, 0, mParticalNum);

  // 画牛奶
  glBindFramebuffer(GL_FRAMEBUFFER, 0);
  glClearColor(1.0f, 0.0f, 0.0f, 1.0f);
  glDisable(GL_DEPTH_TEST);
  glClear(GL_COLOR_BUFFER_BIT);

  mMilkShader->Use();
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, mTextureSdf);
  glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

  // 显示FPS
  char title[128] = "";
  sprintf_s(title, "Fluid Simulation FPS=%.2f", CalculateFPS());
  glfwSetWindowTitle(mWindow, title);

  // 确保切换回默认缓冲（如果之前使用了 FBO）
  glBindFramebuffer(GL_FRAMEBUFFER, 0);
  // 恢复视口到全窗口（避免影响后续操作）
  glViewport(0, 0, windowWidth, windowHeight);
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

// 用于捕获鼠标在窗口中的点击等事件，不加会导致点击后窗口无响应
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

void RenderWidget::LoadVertexes_new(std::vector<glm::vec2>& mPositions, std::vector<float>& mDensity) {
  glBindVertexArray(mVaoParticals);
  glBindBuffer(GL_ARRAY_BUFFER, mPositionBuffer);
  glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec2) * mPositions.size(), mPositions.data(), GL_STATIC_DRAW);
  glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
  glEnableVertexAttribArray(0);

  glBindBuffer(GL_ARRAY_BUFFER, mDensityBuffer);
  glBufferData(GL_ARRAY_BUFFER, sizeof(float) * mDensity.size(), mDensity.data(), GL_STATIC_DRAW);
  glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof(float), (void*)0);
  glEnableVertexAttribArray(1);
  glBindVertexArray(0);

  mParticalNum = mPositions.size();
}

bool RenderWidget::CreateRenderWindow() {
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    

    // 创建窗口
    mWindow = glfwCreateWindow(mWindowWidth, mWindowHeight, "Fluid Simulation", NULL, NULL);
    if (mWindow == nullptr)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return false;
    }
    glfwSetWindowPos(mWindow, 100, 100);
    glfwMakeContextCurrent(mWindow);

    // 注册回调函数
    glfwSetWindowUserPointer(mWindow, this);
    glfwSetFramebufferSizeCallback(mWindow,  ResizeCallback);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return false;
    }

    // 初始化 ImGui
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();

    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;  // 启用键盘控制
    // 设置 ImGui 风格（可选）
    ImGui::StyleColorsDark();

    ImGui_ImplGlfw_InitForOpenGL(mWindow, true);
    ImGui_ImplOpenGL3_Init("#version 130");
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
    // 找到this指针
    auto thisPtr = reinterpret_cast<RenderWidget*>(glfwGetWindowUserPointer(window));
    glViewport(0, 0, width, height);
}

void RenderWidget::DrawImGuiSidebar() {
  // 全局变量存储参数（或封装在类中）
  float g_Scale = 1.0f;
  bool g_EnableEffect = false;
  glm::vec3 g_Color(1.0f, 1.0f, 1.0f);
  // 获取窗口尺寸
  int windowWidth, windowHeight;
  glfwGetWindowSize(mWindow, &windowWidth, &windowHeight);
  // 渲染 ImGui 前重置视口到全窗口
  glViewport(0, 0, windowWidth, windowHeight);
  // 设置侧边栏位置和大小（右侧 20% 宽度）
  float sidebarWidth = windowWidth * 0.2f;
  ImGui::SetNextWindowPos(ImVec2(windowWidth - sidebarWidth, 0), ImGuiCond_Always);
  ImGui::SetNextWindowSize(ImVec2(sidebarWidth, windowHeight), ImGuiCond_Always);

  // 开始绘制窗口（禁用滚动条和调整大小）
  ImGui::Begin("Control Panel", nullptr, ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoScrollbar);
  // 添加按钮
  if (ImGui::Button("Reset Parameters")) {
    // 按钮点击事件
    g_Scale = 1.0f;
    g_EnableEffect = false;
    g_Color = glm::vec3(1.0f);
  }
  // 添加滑动条
  ImGui::SliderFloat("Scale", &g_Scale, 0.1f, 2.0f);
  ImGui::Checkbox("Enable Effect", &g_EnableEffect);
  // 颜色选择器
  ImGui::ColorEdit3("Color", &g_Color[0]);
  // 结束窗口
  ImGui::End();
  // 提交 ImGui 数据到 GPU
  ImGui::Render();
  glDisable(GL_DEPTH_TEST);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glViewport(0, 0, windowWidth, windowHeight);
  ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
  glfwSwapBuffers(mWindow);
}