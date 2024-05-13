#include <iostream>
#include <string>
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"

using namespace Pythia8;

int countEventsInLHE(const std::string& filename) {
    std::ifstream file(filename);
    std::string line;
    int count = 0;
    while (getline(file, line)) {
        if (line.find("<event>") != std::string::npos) {
            ++count;
        }
    }
    return count;
}

int main(int argc, char* argv[]) {

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " input.lhe -o output.hepmc" << std::endl;
        return 1;
    }
    
    std::string inputFileName = argv[1];
    std::string outputFileName;

    // Parse command line arguments for output file
    for (int i = 2; i < argc; ++i) {
        if (std::string(argv[i]) == "-o" && i + 1 < argc) {
            outputFileName = argv[i + 1];
            break;
        }
    }

    if (outputFileName.empty()) {
        std::cerr << "Output file not specified." << std::endl;
        return 1;
    }

    int numEvents = countEventsInLHE(inputFileName);
    std::cout << "Number of events in LHE file: " << numEvents << std::endl;
    Pythia pythia;
    // pythia.readFile( inputFileName);

    // 从 LHE 文件读取事件
    pythia.readString("Print:quiet = off");
    pythia.readString("Beams:frameType = 4");
    pythia.readString("Beams:LHEF = " + inputFileName);
    pythia.readString("Check:epTolErr = 1e-2");  // 提高容错范围
    pythia.readString("Init:showAllSettings = on");  // 显示所有设置
    pythia.readString("Next:numberCount = 1000");    // 每处理 1000 个事件打印状态
    pythia.readString("Next:numberShowInfo = 1");    // 显示首个事件的信息
    pythia.readString("Next:numberShowProcess = 1"); // 显示首个事件的过程信息
    pythia.readString("Next:numberShowEvent = 1");   // 显示首个事件的详细信息
    // W boson decay settings
    // pythia.readString("24:onMode = off");  // close all w decay modes 
    // pythia.readString("24:onIfAny = 11 12");  // open decay channel W->e+nu_e
    // pythia.readString("24:onIfAny = 13 14");  // open decay channel W->mu+nu_mu
    // pythia.readString("24:onIfAny = 15 16"); 
    // 开启 parton shower 和 FS
    pythia.readString("PartonLevel:ISR = on");  // 初始状态辐射
    pythia.readString("PartonLevel:FSR = on");  // 最终状态辐射
    // 对于 e+ e- 碰撞，通常不需要开启多部分子相互作用
    pythia.readString("PartonLevel:MPI = off");  // 多部分子相互作用关闭

    // 初始化 Pythia
    pythia.init();
    // int iEvent = 0;

    HepMC::Pythia8ToHepMC toHepMC;
    HepMC::IO_GenEvent ascii_io(outputFileName, std::ios::out);
    // int numEvents = 10000;
    
    for (int i = 0; i < numEvents; ++i) {
        if (!pythia.next()) {
                        // Report the reason for the failure
            std::cerr << "Failed to generate event #" << i << std::endl;
            // std::cerr << "Status: " << pythia.info.status() << std::endl;
            // std::cerr << "Error code: " << pythia.info.errorCode() << std::endl;
            // std::cerr << "Error message: " << pythia.info.message() << std::endl;

            // Optionally, you might want to break or continue depending on your needs
            continue;
        }

        // 创建一个 HepMC 事件
        HepMC::GenEvent* hepmc_evt = new HepMC::GenEvent();
        toHepMC.fill_next_event(pythia, hepmc_evt);

        // 写入 HepMC 文件
        ascii_io << hepmc_evt;
        delete hepmc_evt;
    }


    // 完成模拟
    pythia.stat();
    return 0;
}
