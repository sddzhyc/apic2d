#ifndef PARTICAL_SYSTEM_H
#define PARTICAL_SYSTEM_H

#include <vector>
#include <list>
#include <glm/glm.hpp>
#include "Global.h"
#include "Parameters.h"

namespace Fluid2d {

    struct NeighborInfo {
        int index;
        float distance;
        float distance2;
        glm::vec2 radius;
    };

    class ParticalSystem {
    public:
        ParticalSystem();
        ~ParticalSystem();

        void SetContainerSize(glm::vec2 corner, glm::vec2 size);

        int32_t AddFluidBlock(glm::vec2 corner, glm::vec2 size, glm::vec2 v0, float particalSpace);

        void SearchNeighbors();

        size_t GetBlockIdByPosition(glm::vec2 position);

        void BuildBlockStructure();

    private:
        int32_t AddBoundary(glm::vec2 corner, glm::vec2 size);

    public:
        // ���Ӳ���
        float mSupportRadius = 0.025;    // ֧�Ű뾶
        float mSupportRadius2 = mSupportRadius * mSupportRadius;
        float mParticalRadius = 0.005;   // ���Ӱ뾶
        float mParticalDiameter = 2.0 * mParticalRadius;
        float mVolume = 0.8 * mParticalDiameter * mParticalDiameter;    // ���
        float mMass = Para::density0 * mVolume;  // ����
        float mViscosity = 0.01;            // ճ��ϵ��
        float mExponent = 7.0f;              // ѹ��ָ��
        int mStiffness = 50.0f;            // �ն�

        int mStartIndex = 0;
        std::vector<glm::vec2> mPositions;
        std::vector<glm::vec2> mAccleration;
        std::vector<glm::vec2> mVelocity;
        std::vector<float> mDensity;
        std::vector<float> mPressure;
        std::vector<std::vector<NeighborInfo>> mNeighbors;

        // ��������
        glm::vec2 mLowerBound = glm::vec2(-1.0f, -1.0f);
        glm::vec2 mUpperBound = glm::vec2(1.0f, 1.0f);
        glm::vec2 mContainerCenter = glm::vec2(0.0f, 0.0f);
        std::vector<std::vector<int>> mBlocks;
        glm::vec2 mBlockSize = glm::vec2(0.5, 0.5);
        uint32_t mBlockRowNum = 4;
        uint32_t mBlockColNum = 4;

    };

}



#endif // !PARTICAL_SYSTEM_H

