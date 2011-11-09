.test <- function(dir, pattern = "^test_.*\\.R$") {

    require("RUnit");
    package.name <- "NanoStringNorm";

	print(
		paste(
			".test: ",
			system.file("unitTests", package = package.name),
			sep = ""
			)
		);

    suite <- defineTestSuite(
        name = paste(package.name, "RUnit Tests", sep = " "),
        dirs = system.file("unitTests", package = package.name),
        testFileRegexp = pattern,
        rngKind = "default",
        rngNormalKind = "default"
        );
    
    result <- runTestSuite(suite);
    
    printTextProtocol(result, showDetails=TRUE);

    return(result);
    }
