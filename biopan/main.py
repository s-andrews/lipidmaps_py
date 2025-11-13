import os
import logging

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('biopan.log'),
        logging.StreamHandler()
    ]
)

def main():
    # Basic functionality for CLI
    current_dir = os.path.dirname(os.path.abspath(__file__))
    input_file = os.path.join(current_dir, 'tests', 'inputs', 'biopan_test_input.csv')
    
    try:
        # Initialize DataManager
        from src.data_manager import DataManager
        manager = DataManager()
        # Process file (complete pipeline)
        dataset = manager.process_csv(input_file)
        # Extract and print reaction IDs if present in dataset
        reaction_ids = []
        if hasattr(dataset, 'reactions'):
            reaction_ids = [r['reaction_id'] for r in dataset.reactions]
        print(f"Reaction IDs: {reaction_ids}")

        
    except Exception as e:
        logging.error(f"Processing failed: {e}")
        raise


if __name__ == "__main__":
    main()
