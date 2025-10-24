"""
Clean neural network with two hidden layers.
This file family includes guideline-specific variants (G1, G3, G6, G7) and an ALL version.
"""

import numpy as np

# G7 applied: Use bulk operations (vectorized prediction) to reduce Python overhead.
# Why: Fewer Python calls and better BLAS utilization.

class TwoHiddenLayerNeuralNetwork:
    def __init__(self, input_array: np.ndarray, output_array: np.ndarray) -> None:
        """
        Initialize with random weights.
        input_array: training inputs (n_samples x n_features)
        output_array: expected outputs (n_samples x 1)
        """
        self.input_array = input_array
        self.output_array = output_array

        rng = np.random.default_rng(0)
        n_features = self.input_array.shape[1]
        self.input_layer_and_first_hidden_layer_weights = rng.random((n_features, 4))
        self.first_hidden_layer_and_second_hidden_layer_weights = rng.random((4, 3))
        self.second_hidden_layer_and_output_layer_weights = rng.random((3, 1))

        # Buffers for layers
        self.layer_between_input_and_first_hidden_layer = None
        self.layer_between_first_hidden_layer_and_second_hidden_layer = None
        self.layer_between_second_hidden_layer_and_output = None
        self.output = None

    def feedforward(self) -> np.ndarray:
        """Run a forward pass and return the network output."""
        self.layer_between_input_and_first_hidden_layer = sigmoid(
            np.dot(self.input_array, self.input_layer_and_first_hidden_layer_weights)
        )
        self.layer_between_first_hidden_layer_and_second_hidden_layer = sigmoid(
            np.dot(
                self.layer_between_input_and_first_hidden_layer,
                self.first_hidden_layer_and_second_hidden_layer_weights,
            )
        )
        self.layer_between_second_hidden_layer_and_output = sigmoid(
            np.dot(
                self.layer_between_first_hidden_layer_and_second_hidden_layer,
                self.second_hidden_layer_and_output_layer_weights,
            )
        )
        self.output = self.layer_between_second_hidden_layer_and_output
        return self.output

    def back_propagation(self) -> None:
        """
        Standard backprop using sigmoid derivatives.
        """
        # Output layer error
        error_term = 2.0 * (self.output_array - self.layer_between_second_hidden_layer_and_output) * sigmoid_derivative(
            self.layer_between_second_hidden_layer_and_output
        )
        updated_second_hidden_layer_and_output_layer_weights = np.dot(
            self.layer_between_first_hidden_layer_and_second_hidden_layer.T,
            error_term,
        )

        # Backprop to 2nd hidden
        second_hidden_error = np.dot(
            error_term, self.second_hidden_layer_and_output_layer_weights.T
        ) * sigmoid_derivative(self.layer_between_first_hidden_layer_and_second_hidden_layer)
        updated_first_hidden_layer_and_second_hidden_layer_weights = np.dot(
            self.layer_between_input_and_first_hidden_layer.T,
            second_hidden_error,
        )

        # Backprop to 1st hidden
        first_hidden_error = np.dot(
            second_hidden_error, self.first_hidden_layer_and_second_hidden_layer_weights.T
        ) * sigmoid_derivative(self.layer_between_input_and_first_hidden_layer)
        updated_input_layer_and_first_hidden_layer_weights = np.dot(
            self.input_array.T, first_hidden_error
        )

        # Gradient ascent step (as in original)
        self.input_layer_and_first_hidden_layer_weights += updated_input_layer_and_first_hidden_layer_weights
        self.first_hidden_layer_and_second_hidden_layer_weights += updated_first_hidden_layer_and_second_hidden_layer_weights
        self.second_hidden_layer_and_output_layer_weights += updated_second_hidden_layer_and_output_layer_weights

    def train(self, output: np.ndarray, iterations: int, give_loss: bool) -> None:
        """
        Train for 'iterations' steps. If give_loss is True, print loss per step.
        """
        for iteration in range(1, iterations + 1):
            self.output = self.feedforward()
            self.back_propagation()
            if give_loss:
                loss = np.mean(np.square(output - self.output))
                print(f"Iteration {iteration} Loss: {loss}")

    def predict(self, input_arr: np.ndarray):
        """
        (G7) Vectorized prediction: accept 1D or 2D input and compute in bulk.
        Returns an int for a single sample, or a list of ints for a batch.
        """
        arr = np.asarray(input_arr, dtype=float)
        if arr.ndim == 1:
            arr = arr.reshape(1, -1)
        l1 = sigmoid(np.dot(arr, self.input_layer_and_first_hidden_layer_weights))
        l2 = sigmoid(np.dot(l1, self.first_hidden_layer_and_second_hidden_layer_weights))
        out = sigmoid(np.dot(l2, self.second_hidden_layer_and_output_layer_weights)).ravel()
        preds = (out > 0.6).astype(int)
        return int(preds[0]) if preds.size == 1 else preds.tolist()

def sigmoid(value: np.ndarray) -> np.ndarray:
    """Sigmoid activation function."""
    return 1.0 / (1.0 + np.exp(-value))


def sigmoid_derivative(activated: np.ndarray) -> np.ndarray:
    """Derivative of sigmoid; expects activated values (sigmoid(x))."""
    return activated * (1.0 - activated)

def example() -> int:
    """
    Fixed-input example to avoid interactive I/O.
    """
    test_input = np.array(([0, 0, 0], [0, 1, 0], [0, 0, 1]), dtype=np.float64)
    output = np.array(([0], [1], [1]), dtype=np.float64)

    nn = TwoHiddenLayerNeuralNetwork(test_input, output)
    nn.train(output=output, iterations=10, give_loss=False)
    return nn.predict(np.array(([1, 1, 1]), dtype=np.float64))  # still fixed input; API now supports batches


if __name__ == "__main__":
    example()
