#ifndef EVENTRANDOMGEN_H
#define EVENTRANDOMGEN_H

#include <string>
#include "core/Event.h"

/**
 * @class EventRandomGen
 * @brief 基于多重性和 Parton 分布直方图，随机生成一组 Parton 事件。
 */
class EventRandomGen {
public:
    /// Sampling mode for event generation
    enum SamplingMode {
        kToyMode,
        kSampleFromFile
    };
    /**
     * 构造函数
     * @param histFilePath Parton 直方图 ROOT 文件路径。
     *        如果设置环境变量 CW_COAL_PARTON_HIST 则优先使用该值；
     *        否则默认使用安装目录下 DATA_INSTALL_DIR/dist_parton_afART.root
     */
    explicit EventRandomGen(
        const std::string& histFilePath = std::string(DATA_INSTALL_DIR) + "/dist_parton_afART.root"
    );

    /**
     * 生成一个事件。
     * @param nParts 粒子数；如果小于0，则从多重性直方图随机抽取。
     * @param sumBaryonNumber 期望的总重子数（单位为1，默认为0）。
     * @param mode 采样模式，默认为 kSampleFromFile。
     */
     void GenerateEvent(Event& out,
                        int nParts = -1,
                        int sumBaryonNumber = 0,
                        SamplingMode mode = kSampleFromFile) const;

private:

    std::string m_histFilePath;
};

#endif // EVENTRANDOMGEN_H
