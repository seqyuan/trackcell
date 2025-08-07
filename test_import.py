#!/usr/bin/env python3
"""
Simple test script to verify that the trackcell package can be imported correctly.
"""

try:
    import trackcell
    print(f"✓ TrackCell version {trackcell.__version__} imported successfully!")
    
    # Test importing submodules
    import trackcell.io
    print("✓ trackcell.io module imported successfully!")
    
    # Test importing the main function
    from trackcell.io import read_hd_cellseg
    print("✓ read_hd_cellseg function imported successfully!")
    
    print("\n🎉 All imports successful! The package is ready to use.")
    
except ImportError as e:
    print(f"❌ Import error: {e}")
    exit(1)
except Exception as e:
    print(f"❌ Unexpected error: {e}")
    exit(1) 