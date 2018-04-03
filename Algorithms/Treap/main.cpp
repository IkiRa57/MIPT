#include <bits/stdc++.h>

template<typename T>
class Treap {
public:
private:
    enum ORDER {
        ASC, DESC
    };

    struct Node {
        bool reversed, asc, desc;
        size_t count;
        T value, l_value, r_value;
        T add;
        bool assign_flag;
        T assign_value;
        T sum;
        int priority;
        Node *left, *right, *parent;

        Node() :
                reversed(false),
                count(1),
                asc(true),
                desc(true),
                value(0),
                l_value(0),
                r_value(0),
                add(0),
                assign_flag(false),
                assign_value(0),
                sum(0),
                priority(rand()),
                left(nullptr),
                right(nullptr),
                parent(nullptr) {}

        Node(const T &x) :
                reversed(false),
                asc(true),
                desc(true),
                count(1),
                value(x),
                l_value(x),
                r_value(x),
                add(0),
                assign_flag(false),
                assign_value(0),
                sum(x),
                priority(rand()),
                left(nullptr),
                right(nullptr),
                parent(nullptr) {}
    };

    bool comp(const T &a, const T &b, ORDER order) const {
        return order == ASC ? a < b : a > b;
    }

    int cnt(Node *&root) const {
        return (root) ? root->count : 0;
    }

    T sm(Node *&root) const {
        if (!root) {
            return 0;
        }
        if (root->assign_flag) {
            return root->count * root->assign_value;
        }
        return root->sum + root->add * root->count;
    }

    T val(Node *root) const {
        if (root->assign_flag)
            return root->assign_value;
        return root->value + root->add;
    }

    T get_lr_value(Node *root, char side) const {
        if (root->assign_flag)
            return root->assign_value;
        return ((root->reversed ^ (side == 'l')) ? root->l_value : root->r_value) + root->add;
    }

    T get_r_value(Node *root) const {
        return get_lr_value(root, 'r');
    }

    T get_l_value(Node *root) const {
        return get_lr_value(root, 'l');
    }

    bool is_ordered_2(Node *&root, ORDER order) const {
        if (root->assign_flag)
            return true;
        if (!root->reversed) {
            if (order == DESC) {
                return root->desc;
            }
            return root->asc;
        }
        if (order == DESC) {
            return root->asc;
        }
        return root->desc;
    }

    bool is_ordered(Node *&root, ORDER order) const {
        if (!root)
            return true;
        if (root->left) {
            if (!is_ordered_2(root->left, order))
                return false;
            if (comp(root->value, get_r_value(root->left), order))
                return false;
        }
        if (root->right) {
            if (!is_ordered_2(root->right, order))
                return false;
            if (comp(get_l_value(root->right), root->value, order))
                return false;
        }
        return true;
    }

    bool is_desc(Node *&root) const {
        return is_ordered(root, DESC);
    }

    bool is_asc(Node *&root) const {
        return is_ordered(root, ASC);
    }

    void recalc(Node *&root) {
        if (!root)
            return;
        if (root->add != 0 || root->assign_flag || root->reversed) {
            exit(27);
        }
        root->count = 1 + cnt(root->left) + cnt(root->right);
        root->sum = root->value + sm(root->left) + sm(root->right);
        root->l_value = root->value;
        root->r_value = root->value;
        if (root->left) {
            root->left->parent = root;
            root->l_value = get_l_value(root->left);
        }
        if (root->right) {
            root->right->parent = root;
            root->r_value = get_r_value(root->right);
        }
        root->asc = is_asc(root);
        root->desc = is_desc(root);
    }

    void push_assign(Node *root) {
        root->add = 0;
        if (root->left) {
            root->left->add = 0;
            root->left->assign_value = root->assign_value;
            root->left->assign_flag = root->assign_flag;
        }
        if (root->right) {
            root->right->add = 0;
            root->right->assign_value = root->assign_value;
            root->right->assign_flag = root->assign_flag;
        }
        root->sum = root->count * root->assign_value;
        root->asc = root->desc = 1;
        root->value = root->assign_value;
        root->l_value = root->assign_value;
        root->r_value = root->assign_value;
        root->assign_flag = false;
        root->assign_value = 0;
    }

    void push_add(Node *root) {
        root->sum += root->add * root->count;
        root->value += root->add;
        root->l_value += root->add;
        root->r_value += root->add;
        if (root->left) {
            if (root->left->assign_flag) {
                root->left->assign_value += root->add;
                root->left->add = 0;
            } else {
                root->left->add += root->add;
            }
        }
        if (root->right) {
            if (root->right->assign_flag) {
                root->right->assign_value += root->add;
                root->right->add = 0;
            } else {
                root->right->add += root->add;
            }
        }
        root->add = 0;
    }

    void push_reversed(Node *root) {
        root->reversed = false;
        std::swap(root->left, root->right);
        std::swap(root->asc, root->desc);
        std::swap(root->l_value, root->r_value);
        if (root->left) {
            root->left->reversed ^= true;
        }
        if (root->right) {
            root->right->reversed ^= true;
        }
    }

    void Push(Node *root) {
        if (!root)
            return;
        if (root->assign_flag) {
            push_assign(root);
        }
        if (root->add != 0) {
            push_add(root);
        }
        if (root->reversed) {
            push_reversed(root);
        }
    }

    Node *merge(Node *L, Node *R) {
        Push(L);
        Push(R);
        if (!L || !R) {
            if (L) {
                recalc(L);
                return L;
            } else {
                recalc(R);
                return R;
            }
        }
        if (L->priority > R->priority) {
            L->right = merge(L->right, R);
            if (L->right)
                L->right->parent = L;
            L->parent = nullptr;
            recalc(L);
            return L;
        } else {
            R->left = merge(L, R->left);
            if (R->left)
                R->left->parent = R;
            R->parent = nullptr;
            recalc(R);
            return R;
        }
    }

    void split(Node *root, size_t pos, Node *&L, Node *&R) {
        if (!root) {
            L = R = nullptr;
            return;
        }
        Push(root);
        size_t curr_pos = cnt(root->left);
        if (pos <= curr_pos) {
            split(root->left, pos, L, root->left);
            if (root->left)
                root->left->parent = root;
            if (L)
                L->parent = nullptr;
            R = root;
        } else {
            split(root->right, pos - curr_pos - 1, root->right, R);
            if (root->right)
                root->right->parent = root;
            if (R)
                R->parent = nullptr;
            L = root;
        }
        recalc(root);
    }

    int get_ordered_suf(Node *&root, ORDER order) {
        if (!root)
            return -1;
        Push(root);
        Push(root->left);
        Push(root->right);
        if ((order == ASC && root->asc) || (order == DESC && root->desc))
            return 0;
        if (root->right) {
            bool flag = (order == ASC) ? root->right->asc : root->right->desc;
            if (!flag || comp(root->right->l_value, root->value, order))
                return cnt(root->left) + get_ordered_suf(root->right, order) + 1;
        }
        if (!root->left)
            return 0;
        if (comp(root->value, root->left->r_value, order))
            return cnt(root->left);
        return get_ordered_suf(root->left, order);
    }

    int pos_in_ordered(Node *&root, const T &val, ORDER order) {
        if (!root) {
            return -1;
        }
        Push(root);
        Push(root->left);
        Push(root->right);
        if (!comp(root->value, val, order))
            return pos_in_ordered(root->left, val, order);
        if (root->right && comp(root->right->l_value, val, order))
            return cnt(root->left) + 1 + pos_in_ordered(root->right, val, order);
        return cnt(root->left);
    }

    template<typename Lambda>
    void apply_on_segment(Node *root, size_t l, size_t r, Lambda &&func) {
        Node *L, *M, *R;
        split(root, l, L, R);
        split(R, r - l, M, R);
        func(M);
        R = merge(M, R);
        root = merge(L, R);
    }

    void insert(Node *&root, Node *it, size_t pos) {
        Node *L, *R;
        if (!root) {
            root = it;
            root->count = 1;
            root->parent = nullptr;
            return;
        }
        split(root, pos, L, R);
        L = merge(L, it);
        root = merge(L, R);
    }

    void erase(Node *&root, size_t index) {
        if (!root) {
            return;
        }
        Push(root);
        size_t currIndex = cnt(root->left);
        if (currIndex == index) {
            root = merge(root->left, root->right);
            return;
        }
        if (currIndex > index) {
            erase(root->left, index);
        }
        if (currIndex < index) {
            erase(root->right, index - currIndex - 1);
        }
        recalc(root);
    }

    void reverse(Node *&root, size_t l, size_t r) {
        apply_on_segment(root, l, r, [](Node *&vertex) {
            vertex->reversed ^= 1;
        });
    }

    void add_on_segment(Node *&root, size_t l, size_t r, T val) {
        apply_on_segment(root, l, r, [=](Node *&vertex) {
            vertex->add = val;
        });
    }

    void assign_on_segment(Node *&root, size_t l, size_t r, T val) {
        apply_on_segment(root, l, r, [=](Node *&vertex) {
            vertex->add = 0;
            vertex->assign_flag = true;
            vertex->assign_value = val;
        });
    }

    T sum_on_segment(Node *&root, int l, int r) {
        T res = 0;
        apply_on_segment(root, l, r, [&](Node *&vertex) {
            res = vertex->sum;
        });
        return res;
    }

    void ordered_permutation(Node *&root, int l, int r, ORDER order) {
        apply_on_segment(root, l, r, [=](Node *&vertex) {
            if ((order == DESC && vertex->desc) || (order == ASC && vertex->asc)) {
                vertex->reversed ^= true;
                return;
            }
            int pos = get_ordered_suf(vertex, order);
            Node *suf, *begin, *x;
            split(vertex, pos - 1, begin, suf);
            split(suf, 1, x, suf);
            pos = pos_in_ordered(suf, x->value, order);
            Node *l_suf, *r_suf, *y;
            split(suf, pos, l_suf, r_suf);
            split(r_suf, 1, y, r_suf);
            suf = merge(l_suf, x);
            suf = merge(suf, r_suf);
            suf->reversed ^= 1;
            vertex = merge(y, suf);
            vertex = merge(begin, vertex);

        });
    }

    void next_permutation(Node *&root, int l, int r) {
        return ordered_permutation(root, l, r, DESC);
    }

    void prev_permutation(Node *&root, int l, int r) {
        return ordered_permutation(root, l, r, ASC);
    }

    Node *tree_by_node(Node *root) const {
        Node *a = root;
        while (a && a->parent) {
            a = a->parent;
        }
        return a;
    }

    int ind_by_node(Node *root, Node *node) {
        Node *a = node;
        std::stack < Node * > st;
        st.push(a);
        while (a->parent) {
            a = a->parent;
            st.push(a);
        }
        while (!st.empty()) {
            Push(st.top());
            recalc(st.top());
            st.pop();
        }
        a = node;
        int k = cnt(a->left);
        while (a->parent) {
            Node *b = a->parent;
            if (a == b->left) {}
            if (a == b->right) {
                k += (cnt(b->left) + 1);
            }
            a = b;
        }
        return k;
    }

    Node *find_by_ind(Node *&root, int index) const {
        if (index == cnt(root->left)) {
            return root;
        }
        if (index < cnt(root->left)) {
            return find_by_ind(root->left, index);
        } else {
            return find_by_ind(root->right, index - cnt(root->left) - 1);
        }
    }

    void inorderTraversal(Node *root) {
        if (!root) {
            return;
        }
        Push(root);
        inorderTraversal(root->left);
        std::cout << root->value << ' ';
        inorderTraversal(root->right);
    }

    void clean(Node *root) {
        if (!root) {
            return;
        }
        clean(root->left);
        clean(root->right);
        delete root;
    }

public:


    Node *root;

    Treap() : root(nullptr) {}

    ~Treap() {
        clean(root);
    }

    Treap(const T &key) {
        Node *a = new Node(key);
        root = a;
    }

    Treap(Node *rot) {
        root = rot;
    }

    void segment_to_begin(int l, int r) {
        segment_to_begin(root, l, r);
    }

    void reverse(int l, int r) {
        reverse(root, l, r);
    }

    void add_on_segment(int l, int r, T val) {
        add_on_segment(root, l, r, val);
    }

    void assign_on_segment(int l, int r, T val) {
        assign_on_segment(root, l, r, val);
    }

    T sum_on_segment(int l, int r) {
        return sum_on_segment(root, l, r);
    }

    void next_permutation(int l, int r) {
        return next_permutation(root, l, r);
    }

    void prev_permutation(int l, int r) {
        return prev_permutation(root, l, r);
    }

    void insert(const T &key, int pos) {
        Node *x = new Node(key);
        insert(root, x, pos);
    }

    void erase(int pos) {
        erase(root, pos);
    }

    void inorderTraversal() {
        inorderTraversal(root);
    }

    T operator[](int number) {
        return find_by_ind(root, number);
    }
};

int main() {
    std::ios_base::sync_with_stdio(0);
    std::cin.tie(NULL);
    int n, q;
    Treap<long long> a;
    std::cin >> n;
    int x;
    for (int i = 0; i < n; ++i) {
        std::cin >> x;
        a.insert(x, i);
    }
    std::cin >> q;
    int t, l, r, pos;
    for (int i = 0; i < q; ++i) {
        std::cin >> t;
        if (t == 1) {
            std::cin >> l >> r;
            std::cout << a.sum_on_segment(l, r + 1) << std::endl;
        }
        if (t == 2) {
            std::cin >> x >> pos;
            a.insert(x, pos);
        }
        if (t == 3) {
            std::cin >> pos;
            a.erase(pos);
        }
        if (t == 4) {
            std::cin >> x >> l >> r;
            a.assign_on_segment(l, r + 1, x);
        }
        if (t == 5) {
            std::cin >> x >> l >> r;
            a.add_on_segment(l, r + 1, x);
        }
        if (t == 6) {
            std::cin >> l >> r;
            a.next_permutation(l, r + 1);
        }
        if (t == 7) {
            std::cin >> l >> r;
            a.prev_permutation(l, r + 1);
        }
    }
    a.inorderTraversal();
}
