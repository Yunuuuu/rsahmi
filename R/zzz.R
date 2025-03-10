.onAttach <- function(libname, pkgname) {
    version <- packageDescription(pkgname, fields = "Version")
    msg <- paste0(
        "========================================\n",
        pkgname, " version ", version, "\n\n",
        "If you use it in published research, please cite:\n",
        "- Song, Yuxuan PhD; Peng, Yun PhD; Qin, Caipeng PhD; Jiang, Shan PhD; Lin, Jiaxing PhD; Lai, Shicong MD; Wu, Jilin PhD; Ding, Mengting PhD; Du, Yiqing PhD*; Yu, Luping MD*; Xu, Tao MD*. Antibiotic use attenuates response to immune checkpoint blockade in urothelial carcinoma via inhibiting CD74-MIF/COPA: revealing cross-talk between anti-bacterial immunity and ant-itumor immunity. International Journal of Surgery 111(1):p 972-987, January 2025. | DOI: 10.1097/JS9.0000000000001901. \n",
        "- Ghaddar B, Blaser MJ, De S. Denoising sparse microbial signals from single-cell sequencing of mammalian host tissues. Nat Comput Sci. 2023 Sep;3(9):741-747. doi: 10.1038/s43588-023-00507-1. Epub 2023 Sep 18. PMID: 37946872; PMCID: PMC10634611.\n",
        "\n",
        "This message can be suppressed by: ",
        sprintf("suppressPackageStartupMessages(library(%s))\n", pkgname),
        "========================================"
    )
    packageStartupMessage(msg)
}
