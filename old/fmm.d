import std.algorithm;
import std.array;
import std.container.binaryheap;
import std.conv;
import std.format;
import std.math;
import std.stdio;
import std.typecons;

enum State {known, trial, far};

struct Node {
	double value = double.infinity;
	State state = State.far;
	int i;
	int j;

	bool isKnown() const pure {
		return state == State.known;
	}

	bool isTrial() const pure {
		return state == State.trial;
	}

	bool isFar() const pure {
		return state == State.far;
	}

	bool isActive() const pure {
		return state == State.trial || state == State.far;
	}

	string getStateString() const {
		return isKnown() ? "known" : isTrial() ? "trial" : "far";
	}

	string toString() const {
		return format!"Node(%g, %s, %d, %d)"(value, getStateString(), i, j);
	}
}
	
struct Mesh {
	int width;
	int height;
	double h;
	Node[] nodes;

	this(int width, int height, double h) {
		this.width = width;
		this.height = height;
		this.h = h;
		nodes.length = width*height;
		for (int i = 0; i < height; ++i) {
			for (int j = 0; j < width; ++j) {
				Node node;
				node.i = i;
				node.j = j;
				this[i, j] = node;
			}
		}
	}

	string toString() {
		auto app = appender!string();
		for (int i = 0; i < height; ++i) {
			for (int j = 0; j < width - 1; ++j) {
				app.put(to!string(this[i, j]));
				app.put(" ");
			}
			app.put(to!string(this[i, width - 1]));
			app.put('\n');
		}
		return app.data();
	}
	
	ref Node opIndex(int i, int j) {
		return nodes[width*i + j];
	}

	auto getNodesByState(State state) {
		Node*[] tmp;
		foreach (ref Node node; nodes)
			if (node.state == state)
				tmp ~= &node;
		return tmp;
	}

	auto getKnownNodes() {
		return getNodesByState(State.known);
	}

	auto getTrialNodes() {
		return getNodesByState(State.trial);
	}

	auto getActiveNodes() {
		Node*[] tmp;
		foreach (ref Node node; nodes)
			if (node.isActive())
				tmp ~= &node;
		return tmp;
	}

	auto getNeighbors(Node* node) {
		auto const i = node.i;
		auto const j = node.j;
		Node*[] neighbors;
		neighbors.reserve(4);
		if (i - 1 >= 0) neighbors ~= &this[i - 1, j];
		if (i + 1 < height) neighbors ~= &this[i + 1, j];
		if (j - 1 >= 0) neighbors ~= &this[i, j - 1];
		if (j + 1 < width) neighbors ~= &this[i, j + 1];
		return neighbors;
	}

	enum Direction { up, down, left, right };
	
	Node*[Direction] getNeighborsByDirection(Node* node) {
		auto const i = node.i;
		auto const j = node.j;
		Node*[Direction] neighbors;
		if (i - 1 >= 0) neighbors[Direction.up] = &this[i - 1, j];
		if (i + 1 < height) neighbors[Direction.down] = &this[i + 1, j];
		if (j - 1 >= 0) neighbors[Direction.left] = &this[i, j - 1];
		if (j + 1 < width) neighbors[Direction.right] = &this[i, j + 1];
		return neighbors;
	}

	void updateValue(Node* node) {
		Node*[Direction] known;
		foreach (Direction direction, Node* neighbor;
				 getNeighborsByDirection(node)) {
			if (neighbor.isKnown()) {
				known[direction] = neighbor;
			}
		}
		
		Node* n;
		double T, T1, T2, f = 1, rhs = h*h/(f*f);

		if (Direction.up in known || Direction.down in known) {
			if (Direction.left in known || Direction.right in known) {
				if (Direction.up in known) {
					n = known[Direction.up];
					if (n.isKnown()) T1 = n.value;
				}
				if (Direction.down in known) {
					n = known[Direction.down];
					if (n.isKnown()) T1 = min(T1, n.value);
				}

				if (Direction.left in known) {
					n = known[Direction.left];
					if (n.isKnown()) T2 = n.value;
				}
				if (Direction.right in known) {
					n = known[Direction.right];
					if (n.isKnown()) T2 = min(T2, n.value);
				}

				T = (T1 + T2 + sqrt(2*rhs - pow(T1 - T2, 2)))/2;
			} else {
				if (Direction.up in known) {
					n = known[Direction.up];
					if (n.isKnown()) T1 = n.value;
				}
				if (Direction.down in known) {
					n = known[Direction.down];
					if (n.isKnown()) T1 = min(T1, n.value);
				}

				T = rhs + T1;
			}
		} else {
			assert(Direction.left in known || Direction.right in known);

			if (Direction.left in known) {
				n = known[Direction.left];
				if (n.isKnown()) T1 = n.value;
			}
			if (Direction.right in known) {
				n = known[Direction.right];
				if (n.isKnown()) T1 = min(T1, n.value);
			}

			T = rhs + T1;
		}

		node.value = T;
	}

	bool finished() const {
		return all!(node => node.isKnown())(nodes);
	}

	void addBoundaryPoint(int i, int j) {
		this[i, j].value = 0;
		this[i, j].state = State.known;
	}
	
	string toMatlabMatrixString() {
		assert(finished());
		auto app = appender!string();
		app.put('[');
		for (int i = 0; i < height - 1; ++i) {
			for (int j = 0; j < width - 1; ++j) {
				app.put(to!string(this[i, j].value));
				app.put(' ');
			}
			app.put(to!string(this[i, width - 1].value));
			app.put("; ");
		}
		for (int j = 0; j < width - 1; ++j) {
			app.put(to!string(this[height - 1, j].value));
			app.put(' ');
		}
		app.put(to!string(this[height - 1, width - 1].value));
		app.put(']');
		return app.data;
	}
}

void main() {
	Mesh mesh = Mesh(21, 21, 1);
	mesh.addBoundaryPoint(3, 3);
	mesh.addBoundaryPoint(10, 20);
	mesh.addBoundaryPoint(1, 17);
	mesh.addBoundaryPoint(19, 1);

	foreach (node; mesh.getKnownNodes()) {
		foreach (neighbor; mesh.getNeighbors(node)) {
			if (neighbor.isActive()) {
				neighbor.state = State.trial;
				mesh.updateValue(neighbor);
			}
		}
	}

	while (!mesh.finished()) {
		auto node = minElement!"a.value"(mesh.getActiveNodes());
		node.state = State.known;

		foreach (neighbor; mesh.getNeighbors(node)) {
			if (neighbor.isActive()) {
				mesh.updateValue(neighbor);
				neighbor.state = State.trial;
			}
		}
	}

	writeln(mesh.toMatlabMatrixString());
}
