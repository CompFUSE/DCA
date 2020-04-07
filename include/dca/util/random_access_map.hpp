// Copyright (C) 2020 ETH Zurich
// Copyright (C) 2020 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Provides a map with O(log n) insertion, removal and random access (access of the value relative
// to the i-th lowest key. Useful for a random selection in an ordered list with variable size.
// Implemented as an augmented red-black tree.
// Code to fix RB color violation from https://www.geeksforgeeks.org/red-black-tree-set-1-introduction-2/

#ifndef DCA_UTIL_TYPE_RANDOM_ACCESS_MAP_HPP
#define DCA_UTIL_TYPE_RANDOM_ACCESS_MAP_HPP

#include <cassert>
#include <initializer_list>
#include <functional>
#include <stack>
#include <stdexcept>

#include "dca/util/fixed_size_allocator.hpp"

namespace dca {
namespace util {

// Precondition: elements of type Key have full order.
template <class Key, class Value, std::size_t chunk_size = 256>
class RandomAccessMap {
public:
  RandomAccessMap() = default;
  RandomAccessMap(const std::initializer_list<std::pair<Key, Value>>& list);
  RandomAccessMap(const RandomAccessMap& rhs);
  RandomAccessMap(RandomAccessMap&& rhs);
  RandomAccessMap& operator=(const RandomAccessMap& rhs);
  RandomAccessMap& operator=(RandomAccessMap&& rhs);

  ~RandomAccessMap();

  // Insert new key, value pair.
  // Precondition: key is not already in the map.
  void insert(const Key& key, const Value& value);

  // Remove the node relative to key.
  // Precondition: key is in the map.
  void erase(const Key& key);

  // Returns a reference to the value associated with key.
  // Precondition: key is in the map.
  Value& find(const Key& key);
  const Value& find(const Key& key) const;

  // Returns a reference to the value relative to the i-th key.
  // Precondition: 0 <= index < size()
  const Value& operator[](const std::size_t index) const;

  // Number of keys stored in the map.
  std::size_t size() const {
    return root_ ? root_->subtree_size : 0;
  }

  // Returns an array of ordered keys and value pairs.
  std::vector<std::pair<Key, Value>> linearize() const;

  bool checkConsistency() const;

private:
  enum Color : std::uint8_t { RED, BLACK };

  struct Node {
    Node(const Key& k, const Value& v, Node* p) : parent(p), key(k), val(v) {}

    Node* left = nullptr;
    Node* right = nullptr;
    Node* parent = nullptr;

    std::size_t subtree_size = 1;

    Key key;
    Value val;

    Color color = RED;
  };

  void rightRotate(Node* node);

  void leftRotate(Node* node);

  void fixRedRed(Node* node);

  void fixDoubleBlack(Node* node);

  bool isLeftChild(const Node* node) const {
    return node->parent && node->parent->left == node;
  }

  bool isRightChild(const Node* node) const {
    return node->parent && node->parent->right == node;
  }

  auto getUncle(const Node* node) const -> Node*;

  auto getSibling(const Node* node) const -> Node*;

  void updateSubtreeSize(Node* node);

  // moves node down and moves given node in its place
  void moveDown(Node* node, Node* new_parent);

  // Members
  Node* root_ = nullptr;
  dca::util::details::FixedSizeAllocator<Node, chunk_size> allocator_;
};

template <class Key, class Value, std::size_t chunk_size>
RandomAccessMap<Key, Value, chunk_size>::RandomAccessMap(const std::initializer_list<std::pair<Key, Value>>& list) {
  for (const auto& [key, val] : list)
    insert(key, val);
}

template <class Key, class Value, std::size_t chunk_size>
RandomAccessMap<Key, Value, chunk_size>::~RandomAccessMap() {
  std::stack<Node*> to_delete;
  if (root_)
    to_delete.push(root_);

  while (!to_delete.empty()) {
    auto node = to_delete.top();
    to_delete.pop();

    if (node->left)
      to_delete.push(node->left);
    if (node->right)
      to_delete.push(node->right);

    allocator_.destroy(node);
  }
}

template <class Key, class Value, std::size_t chunk_size>
RandomAccessMap<Key, Value, chunk_size>::RandomAccessMap(const RandomAccessMap& rhs) {
  (*this) = rhs;
}

template <class Key, class Value, std::size_t chunk_size>
RandomAccessMap<Key, Value, chunk_size>::RandomAccessMap(RandomAccessMap&& rhs) {
  (*this) = std::move(rhs);
}

template <class Key, class Value, std::size_t chunk_size>
RandomAccessMap<Key, Value, chunk_size>& RandomAccessMap<Key, Value, chunk_size>::operator=(
    const RandomAccessMap<Key, Value, chunk_size>& rhs) {
  if (this != &rhs) {
    *this = std::move(RandomAccessMap());  // clear content.

    const auto linearized = rhs.linearize();
    for (const auto& [key, val] : linearized) {
      insert(key, val);
    }
  }
  return *this;
}

template <class Key, class Value, std::size_t chunk_size>
RandomAccessMap<Key, Value, chunk_size>& RandomAccessMap<Key, Value, chunk_size>::operator=(RandomAccessMap<Key, Value, chunk_size>&& rhs) {
  std::swap(root_, rhs.root_);
  std::swap(allocator_, rhs.allocator_);
  return *this;
}

template <class Key, class Value, std::size_t chunk_size>
void RandomAccessMap<Key, Value, chunk_size>::insert(const Key& key, const Value& val) {
  if (!root_) {
    root_ = allocator_.create(key, val, nullptr);
    root_->color = BLACK;
    return;
  }

  Node* node = root_;
  bool done = false;

  while (!done) {
    if (key == node->key)
      throw(std::logic_error("Key already present."));

    ++node->subtree_size;

    if (key < node->key) {
      if (node->left == nullptr) {
        node->left = allocator_.create(key, val, node);
        done = true;
      }
      node = node->left;
    }
    else {
      if (node->right == nullptr) {
        node->right = allocator_.create(key, val, node);
        done = true;
      }
      node = node->right;
    }
  }

  // Check colors
  fixRedRed(node);

  //  assert(checkConsistency());
}

template <class Key, class Value, std::size_t chunk_size>
void RandomAccessMap<Key, Value, chunk_size>::erase(const Key& key) {
  // Find node to delete
  Node* to_delete = root_;
  while (true) {
    if (!to_delete) {
      throw(std::logic_error("Key not found."));
    }

    if (key == to_delete->key)
      break;

    --to_delete->subtree_size;

    if (key < to_delete->key) {
      to_delete = to_delete->left;
    }
    else {
      to_delete = to_delete->right;
    }
  }

  // to_delete has two children.
  if (to_delete->left != nullptr && to_delete->right != nullptr) {
    Node* const original = to_delete;
    --to_delete->subtree_size;
    to_delete = to_delete->right;
    while (to_delete->left) {
      --to_delete->subtree_size;
      to_delete = to_delete->left;
    }

    // Move key and value from the next in-order node.
    original->key = std::move(to_delete->key);
    original->val = std::move(to_delete->val);
  }

  Node* replacement = to_delete->left ? to_delete->left : to_delete->right;

  auto color = [](const Node* n) { return n ? n->color : BLACK; };
  const bool both_black = color(replacement) == BLACK && to_delete->color == BLACK;
  --to_delete->subtree_size;

  if (both_black) {
    fixDoubleBlack(to_delete);
  }
  else {
    auto sibling = getSibling(to_delete);
    if (sibling && !replacement)
      sibling->color = RED;
    else if (replacement)
      replacement->color = BLACK;
  }

  // delete to_delete from the tree
  Node* parent = to_delete->parent;
  if (isLeftChild(to_delete)) {
    parent->left = replacement;
  }
  else if (isRightChild(to_delete)) {
    parent->right = replacement;
  }
  if (replacement)
    replacement->parent = parent;

  // Update root if necessary.
  if (to_delete == root_) {
    root_ = replacement;
  }

  allocator_.destroy(to_delete);

  //  assert(checkConsistency());
}

template <class Key, class Value, std::size_t chunk_size>
void RandomAccessMap<Key, Value, chunk_size>::fixRedRed(Node* x) {
  // if x is root color it black and return
  if (x == root_) {
    x->color = BLACK;
    return;
  }

  // initialize relatives.
  Node* parent = x->parent;
  Node* grandparent = parent->parent;
  Node* uncle = getUncle(x);

  if (parent->color != BLACK) {
    if (uncle && uncle->color == RED) {
      // uncle is red, perform recoloring and recurse
      parent->color = BLACK;
      uncle->color = BLACK;
      grandparent->color = RED;
      fixRedRed(grandparent);
    }
    else {
      if (isLeftChild(parent)) {
        if (isLeftChild(x)) {
          // for left right
          std::swap(parent->color, grandparent->color);
        }
        else {
          leftRotate(parent);
          std::swap(x->color, grandparent->color);
        }
        // for left left and left right
        rightRotate(grandparent);
      }
      else {
        if (isLeftChild(x)) {
          // for right left
          rightRotate(parent);
          std::swap(x->color, grandparent->color);
        }
        else {
          std::swap(parent->color, grandparent->color);
        }

        // for right right and right left
        leftRotate(grandparent);
      }
    }
  }
}

template <class Key, class Value, std::size_t chunk_size>
const Value& RandomAccessMap<Key, Value, chunk_size>::operator[](const std::size_t index) const {
  if (index >= size())
    throw(std::out_of_range("Index out of range"));

  const Node* node = root_;

  std::size_t on_the_left = 0;
  while (true) {
    assert(node);

    auto new_on_the_left = on_the_left;
    if (node->left)
      new_on_the_left += node->left->subtree_size;

    if (new_on_the_left == index) {
      return node->val;
    }
    else if (new_on_the_left > index) {  // go left
      node = node->left;
    }
    else {  // go right
      on_the_left = new_on_the_left + 1;
      node = node->right;
    }
  }
}

template <class Key, class Value, std::size_t chunk_size>
Value& RandomAccessMap<Key, Value, chunk_size>::find(const Key& key) {
  Node* node = root_;
  while (node) {
    if (node->key == key)
      return node->val;
    else if (key < node->key)
      node = node->left;
    else
      node = node->right;
  }

  throw(std::logic_error("Key not found."));
}

template <class Key, class Value, std::size_t chunk_size>
const Value& RandomAccessMap<Key, Value, chunk_size>::find(const Key& key) const {
  const Node* node = root_;
  while (node) {
    if (node->key == key)
      return node->val;
    else if (key < node->key)
      node = node->left;
    else
      node = node->right;
  }

  throw(std::logic_error("Key not found."));
}

template <class Key, class Value, std::size_t chunk_size>
void RandomAccessMap<Key, Value, chunk_size>::fixDoubleBlack(Node* x) {
  if (x == root_) {  // Reached root
    return;
  }

  Node* sibling = getSibling(x);
  Node* parent = x->parent;

  auto has_red_child = [](Node* n) {
    return (n->left != NULL && n->left->color == RED) || (n->right != NULL && n->right->color == RED);
  };

  if (sibling == NULL) {
    // No sibiling, double black pushed up
    fixDoubleBlack(parent);
  }
  else {
    if (sibling->color == RED) {
      // Sibling red
      parent->color = RED;
      sibling->color = BLACK;
      if (isLeftChild(sibling)) {
        // left case
        rightRotate(parent);
      }
      else {
        // right case
        leftRotate(parent);
      }
      fixDoubleBlack(x);
    }
    else {
      // Sibling black
      if (has_red_child(sibling)) {
        // at least 1 red children
        if (sibling->left != NULL and sibling->left->color == RED) {
          if (isLeftChild(sibling)) {
            // left left
            sibling->left->color = sibling->color;
            sibling->color = parent->color;
            rightRotate(parent);
          }
          else {
            // right left
            sibling->left->color = parent->color;
            rightRotate(sibling);
            leftRotate(parent);
          }
        }
        else {
          if (isLeftChild(sibling)) {
            // left right
            sibling->right->color = parent->color;
            leftRotate(sibling);
            rightRotate(parent);
          }
          else {
            // right right
            sibling->right->color = sibling->color;
            sibling->color = parent->color;
            leftRotate(parent);
          }
        }
        parent->color = BLACK;
      }
      else {
        // 2 black children
        sibling->color = RED;
        if (parent->color == BLACK)
          fixDoubleBlack(parent);
        else
          parent->color = BLACK;
      }
    }
  }
}

template <class Key, class Value, std::size_t chunk_size>
void RandomAccessMap<Key, Value, chunk_size>::rightRotate(Node* const node) {
  // new parent will be node's left child
  Node* new_parent = node->left;

  // update root if current node is root
  if (node == root_)
    root_ = new_parent;

  moveDown(node, new_parent);

  // connect node with new parent's right element
  node->left = new_parent->right;
  // connect new parent's right element with node
  // if it is not nullptr
  if (new_parent->right != nullptr)
    new_parent->right->parent = node;

  // connect new parent with node
  new_parent->right = node;

  updateSubtreeSize(node);
  updateSubtreeSize(new_parent);
}

template <class Key, class Value, std::size_t chunk_size>
void RandomAccessMap<Key, Value, chunk_size>::leftRotate(Node* node) {
  // new parent will be node's right child
  Node* new_parent = node->right;

  // update root_ if current node is root_
  if (node == root_)
    root_ = new_parent;

  moveDown(node, new_parent);

  // connect node with new parent's left element
  node->right = new_parent->left;
  // connect new parent's left element with node
  // if it is not nullptr
  if (new_parent->left != nullptr)
    new_parent->left->parent = node;

  // connect new parent with node
  new_parent->left = node;

  updateSubtreeSize(node);
  updateSubtreeSize(new_parent);
}

template <class Key, class Value, std::size_t chunk_size>
std::vector<std::pair<Key, Value>> RandomAccessMap<Key, Value, chunk_size>::linearize() const {
  std::vector<std::pair<Key, Value>> result;

  std::function<void(const Node*)> helper_func = [&](const Node* node) {
    if (!node)
      return;

    helper_func(node->left);
    result.emplace_back(node->key, node->val);
    helper_func(node->right);
  };

  helper_func(root_);
  assert(result.size() == size());

  return result;
}

template <class Key, class Value, std::size_t chunk_size>
bool RandomAccessMap<Key, Value, chunk_size>::checkConsistency() const {
  bool child_parent_violation = false;
  bool red_red_violation = false;
  bool black_count_violation = false;
  bool subtree_size_violation = false;

  // Returns size of subtree.
  std::function<std::size_t(const Node*)> subtree_size = [&](const Node* node) -> std::size_t {
    if (!node)
      return 0;
    return 1 + subtree_size(node->left) + subtree_size(node->right);
  };

  // Check node consistency and returns number of black nodes in [node, leaves].
  std::function<int(const Node*)> check = [&](const Node* node) {
    if (node == nullptr)
      return 1;

    // Check parent-child relationship.
    if (node->left && node->left->parent != node)
      child_parent_violation = true;
    if (node->right && node->right->parent != node)
      child_parent_violation = true;

    // Check subtree size
    if (node->subtree_size != subtree_size(node))
      subtree_size_violation = true;

    // Check double red
    auto color = [&](const Node* n) { return n ? n->color : BLACK; };
    if (node->color == RED && (color(node->left) == RED || color(node->right) == RED))
      red_red_violation = true;

    // Check black count
    int count_left = check(node->left);
    int count_right = check(node->right);

    if (count_left != count_right)
      black_count_violation = true;

    return count_left + node->color == BLACK ? 1 : 0;
  };

  check(root_);

  return !black_count_violation && !red_red_violation && !child_parent_violation &&
         !subtree_size_violation;
}

template <class Key, class Value, std::size_t chunk_size>
auto RandomAccessMap<Key, Value, chunk_size>::getUncle(const Node* node) const -> Node* {
  const Node* parent = node->parent;
  if (isLeftChild(parent))
    return parent->parent->right;
  else if (isRightChild(parent))
    return parent->parent->left;
  else
    return nullptr;
}

template <class Key, class Value, std::size_t chunk_size>
auto RandomAccessMap<Key, Value, chunk_size>::getSibling(const Node* node) const -> Node* {
  if (isLeftChild(node))
    return node->parent->right;
  else if (isRightChild(node))
    return node->parent->left;
  else
    return nullptr;
}

template <class Key, class Value, std::size_t chunk_size>
void RandomAccessMap<Key, Value, chunk_size>::moveDown(Node* node, Node* new_parent) {
  auto& parent = node->parent;
  if (isLeftChild(node)) {
    parent->left = new_parent;
  }
  else if (isRightChild(node)) {
    parent->right = new_parent;
  }
  new_parent->parent = parent;
  parent = new_parent;
}

template <class Key, class Value, std::size_t chunk_size>
void RandomAccessMap<Key, Value, chunk_size>::updateSubtreeSize(Node* node) {
  node->subtree_size = 1;
  if (node->left)
    node->subtree_size += node->left->subtree_size;
  if (node->right)
    node->subtree_size += node->right->subtree_size;
}

}  // namespace util
}  // namespace dca

#endif  // #define DCA_UTIL_TYPE_RANDOM_ACCESS_MAP_HPP
