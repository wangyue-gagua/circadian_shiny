p1 <- salmonPlotRepCirca("Ghir_D13G018480.1", "Proline-rich spliceosome-associated (PSP) family protein / zinc knuckle (CCHC-type) family protein")
ggsave(p1, "figure/transcriptionPlot/Ghir_D13G018480.1_PSP.pdf", width = 10, height = 10)

p2 <- salmonPlotRepCirca("Ghir_A08G008200.1", "spliceosomal protein U1A")
p3 <- salmonPlotRepCirca("Ghir_D13G018480.2", "Proline-rich spliceosome-associated (PSP) family protein / zinc knuckle (CCHC-type) family protein")

p4 <- salmonPlotRepCirca("Ghir_D06G014080.1", "proline-rich spliceosome-associated (PSP) family protein")
p5 <- salmonPlotRepCirca("Ghir_A13G017740.1", "proline-rich spliceosome-associated (PSP) family protein")

p1 + p2 + p3 + p4 + p5 + plot_layout(ncol = 2, nrow = 3, guides = "collect")
ggsave("figure/transcriptionPlot/PSP_5.pdf", width = 10, height = 10)
