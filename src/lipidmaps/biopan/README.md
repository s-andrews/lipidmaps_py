
# BioPAN Tool

## Overview
The BioPAN framework is designed for managing and processing lipid species data. It provides functionalities to add, manage, and analyze lipid data efficiently.

## Installation
To install the BioPAN framework, clone the repository and install the required dependencies:

```bash
git clone <repository-url>
cd biopan-framework
pip install -r requirements.txt
```

## Usage
After installation, you can use the framework by importing the necessary classes from the `src` package. Here’s a brief example of how to use the `DataManager` class:

```python
from src.data_manager import DataManager

# Initialize DataManager
data_manager = DataManager()

# Add lipid species
data_manager.add_lipid_species(species)

```

## Contributing


## License


## Initial class diagram 
classDiagram
    class BioPANController {
        +runAnalysis(data)
        -ingestor: DataIngestion
        -db: ReactionDatabase
        -analyzer: PathwayAnalyzer
        -visualizer: Visualizer
    }

    class DataIngestion {
        +loadData(filePath)
        +parseLipidData()
    }

    class ReactionDatabase {
        +loadDatabase()
        +compareReactions(lipids)
    }

    class PathwayAnalyzer {
        +analyze(lipids, reactions)
        +identifyPathways()
    }

    class Visualizer {
        +generateGraph(results)
        +exportImage()
    }

    class LipidSpecies {
        +name: str
        +mass: float
        +structure: str
    }

    class Reaction {
        +id: str
        +substrates: list
        +products: list
    }

    class Pathway {
        +name: str
        +reactions: list
    }

    class AnalysisResult {
        +pathways: list
        +score: float
    }

    BioPANController --> DataIngestion
    BioPANController --> ReactionDatabase
    BioPANController --> PathwayAnalyzer
    BioPANController --> Visualizer
    PathwayAnalyzer --> AnalysisResult
    ReactionDatabase --> Reaction
    Pathway --> Reaction
    PathwayAnalyzer --> Pathway
    DataIngestion --> LipidSpecies
