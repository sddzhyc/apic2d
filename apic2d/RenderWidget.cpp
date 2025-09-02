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
    
    // ��ʼ��shader
    //ע��: APIC2D����Ŀ����Ŀ¼�ǣ�D:\Documents\GitHub\apic2d\build\apic2d
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
  // ��ȡ���ڳߴ�
  int windowWidth, windowHeight;
  glfwGetWindowSize(mWindow, &windowWidth, &windowHeight);
  // ����������������ӿڳߴ磨��������ռ���ڿ�ȵ� 20%��
  int mainViewportWidth = static_cast<int>(windowWidth * 0.8f);
  // ��Ĭ��֡���岢�����������ӿ�
  glBindFramebuffer(GL_FRAMEBUFFER, 0);
  glViewport(0, 0, mainViewportWidth, windowHeight);  // �ؼ��������������ӿ�

  float rad0 = 40.0f;
  glm::vec2 minPos(-rad0, -rad0);
  glm::vec2 maxPos(rad0, rad0);
  glm::vec2 dataSize = maxPos - minPos;
      // �������ű�����ȷ��������Ӧ�ӿڸ߶ȣ�
  float targetHeight = windowHeight;
  float targetWidth = mainViewportWidth;
  float scale = std::min(targetWidth / dataSize.x,  // ˮƽ�������ű���
                         targetHeight / dataSize.y  // ��ֱ�������ű���
  );

  // ����ƽ������������ʾ��
  glm::vec2 offset = -0.5f * (minPos + maxPos) * scale;          // ���Ķ���
  offset += glm::vec2(0.5f * targetWidth, 0.5f * targetHeight);  // �ӿ�����

  // �����ӿں�ͶӰ����
  glViewport(0, 0, mainViewportWidth, windowHeight);
  glm::mat4 projection = glm::ortho(0.0f, static_cast<float>(mainViewportWidth),  // ����
                                    0.0f, static_cast<float>(windowHeight),       // �¡���
                                    -1.0f, 1.0f                                   // ����Զƽ��
  );

  // Ӧ�����ź�ƽ�ƣ�ͨ��ģ�;���
  glm::mat4 model = glm::mat4(1.0f);
  model = glm::translate(model, glm::vec3(offset, 0.0f));
  model = glm::scale(model, glm::vec3(scale, scale, 1.0f));
  // �����󴫵ݸ���ɫ��
  mParticalShader->Use();
  mParticalShader->SetMat4("projection", projection);
  mParticalShader->SetMat4("model", model);
  // ��������Ⱦ
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glDrawArrays(GL_POINTS, 0, mParticalNum);
  // ��������Ⱦ�����ݣ���ռ����� 80% ����
  glClearColor(0.3f, 0.3f, 0.3f, 1.0f);
  glEnable(GL_DEPTH_TEST);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  // ������
  glBindFramebuffer(GL_FRAMEBUFFER, 0);
  glClearColor(0.3f, 0.3f, 0.3f, 1.0f);
  glEnable(GL_DEPTH_TEST);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glBindVertexArray(mVaoParticals);
  mParticalShader->Use();
  glEnable(GL_PROGRAM_POINT_SIZE);
  glDrawArrays(GL_POINTS, 0, mParticalNum);

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

  // ȷ���л���Ĭ�ϻ��壨���֮ǰʹ���� FBO��
  glBindFramebuffer(GL_FRAMEBUFFER, 0);
  // �ָ��ӿڵ�ȫ���ڣ�����Ӱ�����������
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

// ���ڲ�������ڴ����еĵ�����¼������ӻᵼ�µ���󴰿�����Ӧ
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
    glfwSetFramebufferSizeCallback(mWindow,  ResizeCallback);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return false;
    }

    // ��ʼ�� ImGui
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();

    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;  // ���ü��̿���
    // ���� ImGui ��񣨿�ѡ��
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
    // �ҵ�thisָ��
    auto thisPtr = reinterpret_cast<RenderWidget*>(glfwGetWindowUserPointer(window));
    glViewport(0, 0, width, height);
}

void RenderWidget::DrawImGuiSidebar() {
  // ȫ�ֱ����洢���������װ�����У�
  float g_Scale = 1.0f;
  bool g_EnableEffect = false;
  glm::vec3 g_Color(1.0f, 1.0f, 1.0f);
  // ��ȡ���ڳߴ�
  int windowWidth, windowHeight;
  glfwGetWindowSize(mWindow, &windowWidth, &windowHeight);
  // ��Ⱦ ImGui ǰ�����ӿڵ�ȫ����
  glViewport(0, 0, windowWidth, windowHeight);
  // ���ò����λ�úʹ�С���Ҳ� 20% ��ȣ�
  float sidebarWidth = windowWidth * 0.2f;
  ImGui::SetNextWindowPos(ImVec2(windowWidth - sidebarWidth, 0), ImGuiCond_Always);
  ImGui::SetNextWindowSize(ImVec2(sidebarWidth, windowHeight), ImGuiCond_Always);

  // ��ʼ���ƴ��ڣ����ù������͵�����С��
  ImGui::Begin("Control Panel", nullptr, ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoScrollbar);
  // ��Ӱ�ť
  if (ImGui::Button("Reset Parameters")) {
    // ��ť����¼�
    g_Scale = 1.0f;
    g_EnableEffect = false;
    g_Color = glm::vec3(1.0f);
  }
  // ��ӻ�����
  ImGui::SliderFloat("Scale", &g_Scale, 0.1f, 2.0f);
  ImGui::Checkbox("Enable Effect", &g_EnableEffect);
  // ��ɫѡ����
  ImGui::ColorEdit3("Color", &g_Color[0]);
  // ��������
  ImGui::End();
  // �ύ ImGui ���ݵ� GPU
  ImGui::Render();
  glDisable(GL_DEPTH_TEST);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glViewport(0, 0, windowWidth, windowHeight);
  ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
  glfwSwapBuffers(mWindow);
}