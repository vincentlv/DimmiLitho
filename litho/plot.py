def plot(ndarray):
    """plots Transversed image, with origin (0,0) at the lower left corner"""
    import matplotlib.pyplot as plt

    plt.imshow(ndarray.T, origin="lower")
    plt.show()
