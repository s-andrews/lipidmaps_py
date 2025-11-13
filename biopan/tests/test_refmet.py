import sys
import os
from src.models.refmet import RefMet
from src.models.sample import QuantifiedLipid
from types import SimpleNamespace

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
current_dir = os.path.dirname(os.path.abspath(__file__))
biopan_test_file = os.path.join(current_dir, 'inputs', 'biopan_small_demo.csv')


import logging

logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Import the actual RefMet class
def test_refmet_annotation():

    # Create test samples
    test_samples = [
        QuantifiedLipid(
            input_name="PC(16:0/18:1)",
            values={"sample1": 1234.5}
        )
    ]
    temp_objs = [SimpleNamespace(sample_name=q.input_name) for q in test_samples]
    logger.info("Calling RefMet.annotate_samples()...")
    result = RefMet.annotate_samples(temp_objs)
    assert result is not None
    print(f"RefMet result: {result}") 
    assert result.get('PC(16:0/18:1)').get('input_name') == "PC(16:0/18:1)"
    assert result.get('PC(16:0/18:1)').get('standardized_name') == "PC 16:0/18:1"
    assert result.get('PC(16:0/18:1)').get('sub_class') == "PC"
