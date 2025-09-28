"""
Stripe Payment Integration for GenoScope
Handles subscriptions, payments, and usage-based billing
"""

import os
import stripe
from typing import Dict, List, Optional, Any
from datetime import datetime, timedelta
from decimal import Decimal
import logging

logger = logging.getLogger(__name__)

# Initialize Stripe
stripe.api_key = os.getenv("STRIPE_SECRET_KEY", "sk_test_...")

class SubscriptionPlans:
    """Available subscription plans"""
    
    FREE = {
        "id": "free",
        "name": "Free",
        "price": 0,
        "stripe_price_id": None,
        "limits": {
            "analyses_per_month": 5,
            "storage_gb": 1,
            "api_calls": 100,
            "vcf_size_mb": 10,
            "features": ["basic_analysis", "clinvar_search"]
        }
    }
    
    PROFESSIONAL = {
        "id": "professional",
        "name": "Professional",
        "price": 99,
        "stripe_price_id": "price_professional_monthly",
        "limits": {
            "analyses_per_month": 50,
            "storage_gb": 100,
            "api_calls": 10000,
            "vcf_size_mb": 500,
            "features": [
                "basic_analysis", "clinvar_search", "gnomad_integration",
                "cosmic_access", "ml_predictions", "batch_processing",
                "priority_queue", "api_access"
            ]
        }
    }
    
    TEAM = {
        "id": "team",
        "name": "Team",
        "price": 499,
        "stripe_price_id": "price_team_monthly",
        "limits": {
            "analyses_per_month": 500,
            "storage_gb": 1000,
            "api_calls": 100000,
            "vcf_size_mb": 2000,
            "team_members": 10,
            "features": [
                "all_professional_features", "team_collaboration",
                "custom_pipelines", "dedicated_support", "sla",
                "audit_logs", "sso_integration"
            ]
        }
    }
    
    ENTERPRISE = {
        "id": "enterprise",
        "name": "Enterprise",
        "price": "custom",
        "stripe_price_id": None,
        "limits": {
            "analyses_per_month": "unlimited",
            "storage_gb": "unlimited",
            "api_calls": "unlimited",
            "vcf_size_mb": "unlimited",
            "team_members": "unlimited",
            "features": [
                "all_features", "on_premise", "custom_development",
                "dedicated_infrastructure", "24_7_support", "training",
                "compliance_reports", "white_label"
            ]
        }
    }
    
    @classmethod
    def get_all_plans(cls) -> List[Dict]:
        """Get all available plans"""
        return [cls.FREE, cls.PROFESSIONAL, cls.TEAM, cls.ENTERPRISE]
    
    @classmethod
    def get_plan(cls, plan_id: str) -> Optional[Dict]:
        """Get plan by ID"""
        plans_map = {
            "free": cls.FREE,
            "professional": cls.PROFESSIONAL,
            "team": cls.TEAM,
            "enterprise": cls.ENTERPRISE
        }
        return plans_map.get(plan_id)


class StripeManager:
    """Manages Stripe operations"""
    
    def __init__(self):
        """Initialize Stripe manager"""
        self.stripe = stripe
        
    def create_customer(self, email: str, name: str, metadata: Dict = None) -> str:
        """
        Create a Stripe customer
        
        Args:
            email: Customer email
            name: Customer name
            metadata: Additional metadata
            
        Returns:
            Stripe customer ID
        """
        try:
            customer = self.stripe.Customer.create(
                email=email,
                name=name,
                metadata=metadata or {}
            )
            logger.info(f"Created Stripe customer: {customer.id}")
            return customer.id
        except stripe.error.StripeError as e:
            logger.error(f"Failed to create customer: {e}")
            raise
    
    def create_subscription(self, 
                          customer_id: str, 
                          price_id: str,
                          trial_days: int = 14) -> Dict:
        """
        Create a subscription for a customer
        
        Args:
            customer_id: Stripe customer ID
            price_id: Stripe price ID
            trial_days: Trial period in days
            
        Returns:
            Subscription details
        """
        try:
            subscription = self.stripe.Subscription.create(
                customer=customer_id,
                items=[{"price": price_id}],
                trial_period_days=trial_days,
                payment_behavior="default_incomplete",
                expand=["latest_invoice.payment_intent"]
            )
            
            return {
                "subscription_id": subscription.id,
                "status": subscription.status,
                "current_period_end": subscription.current_period_end,
                "trial_end": subscription.trial_end,
                "client_secret": subscription.latest_invoice.payment_intent.client_secret
                    if subscription.latest_invoice else None
            }
            
        except stripe.error.StripeError as e:
            logger.error(f"Failed to create subscription: {e}")
            raise
    
    def cancel_subscription(self, subscription_id: str, immediate: bool = False) -> Dict:
        """
        Cancel a subscription
        
        Args:
            subscription_id: Stripe subscription ID
            immediate: Cancel immediately or at period end
            
        Returns:
            Cancellation details
        """
        try:
            if immediate:
                subscription = self.stripe.Subscription.delete(subscription_id)
            else:
                subscription = self.stripe.Subscription.modify(
                    subscription_id,
                    cancel_at_period_end=True
                )
            
            return {
                "status": "cancelled" if immediate else "scheduled_cancellation",
                "cancelled_at": datetime.now().isoformat(),
                "effective_date": datetime.now().isoformat() if immediate 
                    else datetime.fromtimestamp(subscription.current_period_end).isoformat()
            }
            
        except stripe.error.StripeError as e:
            logger.error(f"Failed to cancel subscription: {e}")
            raise
    
    def update_subscription(self, 
                          subscription_id: str, 
                          new_price_id: str) -> Dict:
        """
        Update subscription to a different plan
        
        Args:
            subscription_id: Current subscription ID
            new_price_id: New price ID
            
        Returns:
            Updated subscription details
        """
        try:
            subscription = self.stripe.Subscription.retrieve(subscription_id)
            
            # Update the subscription item
            self.stripe.Subscription.modify(
                subscription_id,
                cancel_at_period_end=False,
                proration_behavior='always_invoice',
                items=[{
                    'id': subscription['items']['data'][0].id,
                    'price': new_price_id,
                }]
            )
            
            return {
                "status": "updated",
                "new_plan": new_price_id,
                "effective_date": datetime.now().isoformat()
            }
            
        except stripe.error.StripeError as e:
            logger.error(f"Failed to update subscription: {e}")
            raise
    
    def create_payment_intent(self, amount: int, currency: str = "usd") -> Dict:
        """
        Create a one-time payment intent
        
        Args:
            amount: Amount in cents
            currency: Currency code
            
        Returns:
            Payment intent details
        """
        try:
            intent = self.stripe.PaymentIntent.create(
                amount=amount,
                currency=currency,
                automatic_payment_methods={"enabled": True}
            )
            
            return {
                "payment_intent_id": intent.id,
                "client_secret": intent.client_secret,
                "amount": amount,
                "currency": currency
            }
            
        except stripe.error.StripeError as e:
            logger.error(f"Failed to create payment intent: {e}")
            raise
    
    def create_usage_record(self, 
                          subscription_item_id: str,
                          quantity: int,
                          timestamp: Optional[int] = None) -> Dict:
        """
        Record usage for metered billing
        
        Args:
            subscription_item_id: Subscription item ID
            quantity: Usage quantity
            timestamp: Unix timestamp
            
        Returns:
            Usage record details
        """
        try:
            usage_record = self.stripe.SubscriptionItem.create_usage_record(
                subscription_item_id,
                quantity=quantity,
                timestamp=timestamp or int(datetime.now().timestamp())
            )
            
            return {
                "id": usage_record.id,
                "quantity": usage_record.quantity,
                "timestamp": usage_record.timestamp
            }
            
        except stripe.error.StripeError as e:
            logger.error(f"Failed to create usage record: {e}")
            raise
    
    def get_customer_portal_url(self, customer_id: str, return_url: str) -> str:
        """
        Generate customer portal URL for billing management
        
        Args:
            customer_id: Stripe customer ID
            return_url: URL to return to after portal session
            
        Returns:
            Portal session URL
        """
        try:
            session = self.stripe.billing_portal.Session.create(
                customer=customer_id,
                return_url=return_url
            )
            return session.url
            
        except stripe.error.StripeError as e:
            logger.error(f"Failed to create portal session: {e}")
            raise


class UsageTracker:
    """Track and enforce usage limits"""
    
    def __init__(self, redis_client=None):
        """
        Initialize usage tracker
        
        Args:
            redis_client: Redis client for caching
        """
        self.redis = redis_client
        
    def track_analysis(self, user_id: str, analysis_type: str, size_mb: float) -> Dict:
        """
        Track an analysis job
        
        Args:
            user_id: User ID
            analysis_type: Type of analysis
            size_mb: Size in MB
            
        Returns:
            Usage update
        """
        key = f"usage:{user_id}:{datetime.now().strftime('%Y-%m')}"
        
        # Increment counters
        if self.redis:
            self.redis.hincrby(key, "analyses_count", 1)
            self.redis.hincrbyfloat(key, "storage_mb", size_mb)
            self.redis.expire(key, 60 * 60 * 24 * 35)  # Expire after 35 days
            
            # Get current usage
            usage = {
                "analyses_count": int(self.redis.hget(key, "analyses_count") or 0),
                "storage_mb": float(self.redis.hget(key, "storage_mb") or 0)
            }
        else:
            # Fallback to in-memory (not recommended for production)
            usage = {
                "analyses_count": 1,
                "storage_mb": size_mb
            }
        
        # Log the usage
        logger.info(f"Usage tracked for user {user_id}: {usage}")
        
        return usage
    
    def check_limits(self, user_id: str, plan: Dict) -> Dict[str, bool]:
        """
        Check if user is within plan limits
        
        Args:
            user_id: User ID
            plan: User's subscription plan
            
        Returns:
            Dict with limit status
        """
        key = f"usage:{user_id}:{datetime.now().strftime('%Y-%m')}"
        
        if self.redis:
            analyses_count = int(self.redis.hget(key, "analyses_count") or 0)
            storage_mb = float(self.redis.hget(key, "storage_mb") or 0)
            api_calls = int(self.redis.hget(key, "api_calls") or 0)
        else:
            analyses_count = storage_mb = api_calls = 0
        
        limits = plan["limits"]
        
        return {
            "can_analyze": analyses_count < limits["analyses_per_month"],
            "can_store": (storage_mb / 1024) < limits["storage_gb"],
            "can_call_api": api_calls < limits["api_calls"],
            "usage": {
                "analyses": f"{analyses_count}/{limits['analyses_per_month']}",
                "storage": f"{storage_mb/1024:.1f}/{limits['storage_gb']} GB",
                "api_calls": f"{api_calls}/{limits['api_calls']}"
            }
        }
    
    def get_usage_summary(self, user_id: str) -> Dict:
        """
        Get usage summary for current month
        
        Args:
            user_id: User ID
            
        Returns:
            Usage summary
        """
        key = f"usage:{user_id}:{datetime.now().strftime('%Y-%m')}"
        
        if self.redis:
            usage = {
                "month": datetime.now().strftime("%B %Y"),
                "analyses_count": int(self.redis.hget(key, "analyses_count") or 0),
                "storage_mb": float(self.redis.hget(key, "storage_mb") or 0),
                "api_calls": int(self.redis.hget(key, "api_calls") or 0),
                "last_analysis": self.redis.hget(key, "last_analysis"),
            }
        else:
            usage = {
                "month": datetime.now().strftime("%B %Y"),
                "analyses_count": 0,
                "storage_mb": 0,
                "api_calls": 0,
                "last_analysis": None
            }
        
        # Calculate costs for metered usage
        usage["estimated_cost"] = self._calculate_usage_cost(usage)
        
        return usage
    
    def _calculate_usage_cost(self, usage: Dict) -> float:
        """
        Calculate usage-based costs
        
        Args:
            usage: Usage data
            
        Returns:
            Estimated cost in USD
        """
        # Pricing model (example)
        PRICE_PER_ANALYSIS = 2.0  # $2 per analysis over limit
        PRICE_PER_GB = 0.10  # $0.10 per GB over limit
        PRICE_PER_1000_API = 1.0  # $1 per 1000 API calls over limit
        
        # This would check against plan limits
        # For now, return a mock calculation
        cost = 0.0
        
        # Add overage charges
        if usage["analyses_count"] > 50:  # Assuming 50 included
            cost += (usage["analyses_count"] - 50) * PRICE_PER_ANALYSIS
        
        if usage["storage_mb"] > 102400:  # 100 GB included
            cost += ((usage["storage_mb"] - 102400) / 1024) * PRICE_PER_GB
        
        if usage["api_calls"] > 10000:  # 10k included
            cost += ((usage["api_calls"] - 10000) / 1000) * PRICE_PER_1000_API
        
        return round(cost, 2)


class BillingWebhooks:
    """Handle Stripe webhooks"""
    
    @staticmethod
    def handle_webhook(payload: bytes, signature: str, webhook_secret: str) -> Dict:
        """
        Handle incoming Stripe webhook
        
        Args:
            payload: Raw request body
            signature: Stripe signature header
            webhook_secret: Webhook endpoint secret
            
        Returns:
            Processing result
        """
        try:
            event = stripe.Webhook.construct_event(
                payload, signature, webhook_secret
            )
            
            # Handle different event types
            handlers = {
                'customer.subscription.created': BillingWebhooks._handle_subscription_created,
                'customer.subscription.updated': BillingWebhooks._handle_subscription_updated,
                'customer.subscription.deleted': BillingWebhooks._handle_subscription_deleted,
                'invoice.payment_succeeded': BillingWebhooks._handle_payment_succeeded,
                'invoice.payment_failed': BillingWebhooks._handle_payment_failed,
            }
            
            handler = handlers.get(event['type'])
            if handler:
                return handler(event)
            
            return {"status": "unhandled", "type": event['type']}
            
        except stripe.error.SignatureVerificationError as e:
            logger.error(f"Webhook signature verification failed: {e}")
            raise
        except Exception as e:
            logger.error(f"Webhook processing failed: {e}")
            raise
    
    @staticmethod
    def _handle_subscription_created(event: Dict) -> Dict:
        """Handle subscription creation"""
        subscription = event['data']['object']
        logger.info(f"Subscription created: {subscription['id']}")
        
        # Update user's subscription in database
        # Send welcome email
        # Grant access to features
        
        return {
            "status": "processed",
            "action": "subscription_activated",
            "subscription_id": subscription['id']
        }
    
    @staticmethod
    def _handle_subscription_updated(event: Dict) -> Dict:
        """Handle subscription update"""
        subscription = event['data']['object']
        logger.info(f"Subscription updated: {subscription['id']}")
        
        # Update user's plan
        # Adjust feature access
        
        return {
            "status": "processed",
            "action": "subscription_updated",
            "subscription_id": subscription['id']
        }
    
    @staticmethod
    def _handle_subscription_deleted(event: Dict) -> Dict:
        """Handle subscription cancellation"""
        subscription = event['data']['object']
        logger.info(f"Subscription cancelled: {subscription['id']}")
        
        # Revoke access
        # Send cancellation email
        # Offer win-back incentive
        
        return {
            "status": "processed",
            "action": "subscription_cancelled",
            "subscription_id": subscription['id']
        }
    
    @staticmethod
    def _handle_payment_succeeded(event: Dict) -> Dict:
        """Handle successful payment"""
        invoice = event['data']['object']
        logger.info(f"Payment succeeded for invoice: {invoice['id']}")
        
        # Update payment status
        # Send receipt
        
        return {
            "status": "processed",
            "action": "payment_received",
            "invoice_id": invoice['id']
        }
    
    @staticmethod
    def _handle_payment_failed(event: Dict) -> Dict:
        """Handle failed payment"""
        invoice = event['data']['object']
        logger.info(f"Payment failed for invoice: {invoice['id']}")
        
        # Send payment failure notification
        # Retry payment
        # Potentially restrict access
        
        return {
            "status": "processed",
            "action": "payment_failed",
            "invoice_id": invoice['id']
        }
